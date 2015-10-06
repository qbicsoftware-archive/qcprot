import os
import sys
import subprocess
import tempfile
import uuid
import shutil
from datetime import datetime
from os.path import join as pjoin
from os.path import exists as pexists
from xml.etree import ElementTree
import hashlib
import base64
import csv
import glob

configfile: "config.json"
workdir: config['var']

SNAKEDIR = config['src']

try:
    VERSION = subprocess.check_output(
        ['git', 'describe', '--tags', '--always', '--dirty'],
        cwd=SNAKEDIR
    ).decode().strip()
except subprocess.CalledProcessError:
    VERSION = 'unknown'

R_HOME = os.path.join(SNAKEDIR, 'r_scripts')

DATA = config['data']
RESULT = config['result']
LOGS = config['logs']
REF = config['ref']
INI_PATH = config['etc']

def data(path):
    return os.path.join(DATA, path)
    
def result(path):
    return os.path.join(RESULT, path)

if 'params' not in config:
    config['params'] = {}

if 'fasta' not in config['params']:
    fastas = glob.glob(os.path.join(config['ref'], '*.fasta'))
    fastas = [os.path.basename(f) for f in fastas]
    if not fastas:
        raise ValueError("no fasta files supplied.")
    config['params']['fasta'] = fastas

try:
    path = subprocess.check_output(["which", "IDMerger"]).decode()
    OPENMS_BIN = os.path.dirname(path)
except subprocess.CalledProcessError:
    OPENMS_BIN = "/usr/bin"

class OpenMS:
    def __init__(self, bin_path, ini_dir, log_dir):
        self._bin_path = bin_path
        self._tools = os.listdir(bin_path)
        self._ini_path = ini_dir
        if not pexists(log_dir):
            os.mkdir(log_dir)
        self._log_dir = log_dir

    def __getattr__(self, name):
        if name not in self._tools:
            raise ValueError("OpenMS tool not found: %s" % name)

        def wrapper(input, output, logfile=None, extra_args=None, ini=None):
            if output is not None and not pexists(os.path.dirname(str(output))):
                os.mkdir(os.path.dirname(str(output)))
            if extra_args is None:
                extra_args = []
            if logfile is None:
                if output is None:
                    identifier = str(uuid.uuid4())
                else:
                    identifier = os.path.basename(str(output))
                logfile = '{}_{}'.format(name, identifier)

            command = [pjoin(self._bin_path, name)]
            if input is not None:
                command += ['-in'] + [str(input)]
            if output is not None:
                command += ['-out'] + [str(output)]
            if ini is not None:
                if isinstance(ini, list):
                    if not len(ini) == 1:
                        raise ValueError("Invalid params.")
                    ini = ini[0]
                command += ['-ini', str(ini).split(',', 1)[1]]
            command += extra_args

            log_std = pjoin(self._log_dir, logfile + '.out')
            log_err = pjoin(self._log_dir, logfile + '.err')
            with open(log_std, 'w') as out:
                with open(log_err, 'w') as err:
                    subprocess.check_call(command, stdout=out, stderr=err)

        return wrapper

openms = OpenMS(OPENMS_BIN, INI_PATH, LOGS)

# store the content of the ini file, so that snakemake will run
# rules agrain if the parameters inside the file change.
def params(name):
    path = pjoin(INI_PATH, name + '.ini')
    try:
        with open(path, 'rb') as f:
            # TODO replace makes sure that there are no wildcards
            # in the file content. This is not exactly clean ;-)
            return "{},{}".format(hashlib.sha256(f.read()).hexdigest(), path)
    except FileNotFoundError as e:
        raise ValueError("ini file '%s' not found" % path) from e


INPUT_FILES = []
for name in os.listdir(DATA):
    if name.lower().endswith('.mzml'):
        INPUT_FILES.append(os.path.basename(name)[:-5])
        if not name.endswith('.mzML'):
            raise ValueError("Extension mzML is case sensitive")
    else:
        print("Ignoring unknown input file %s" % name)


rule all:
    input:
        result("proteins.idXML"),
        result("proteins.csv"),
        result("peptides.csv")


rule FileFilter:
    input: os.path.join(DATA, "{name}.mzML")
    output: "FileFilter/{name}.mzML"
    run:
        openms.FileFilter(input, output, extra_args=['-sort'])


rule InjectionTime:
    input: os.path.join(DATA, "{name}.mzML")
    output: "InjectionTime/{name}.csv"
    run:
        times = []
        mses = []
        retentions = []

        parser = ElementTree.iterparse(input[0], ("start", "end"))
        _, root = next(parser)

        for event, elem in parser:
            if event == "end" and elem.tag == "{http://psi.hupo.org/ms/mzml}spectrum":
                try:
                    injection_time = elem.find('.//{http://psi.hupo.org/ms/mzml}cvParam[@accession="MS:1000927"]').get("value")
                    times.append(injection_time)
                    ms = elem.find("{http://psi.hupo.org/ms/mzml}cvParam[@accession='MS:1000511']").get('value')
                    mses.append(ms)
                    ret = elem.find('.//{http://psi.hupo.org/ms/mzml}cvParam[@accession="MS:1000016"]').get('value')
                    retentions.append(ret)
                except AttributeError:
                    pass
                root.clear()
        with open(output[0], "w") as f:
            f.write(",".join(["time", "mslevel", "rt"]) + "\n")
            for data in zip(times, mses, retentions):
                f.write(",".join(data) + "\n")


rule CombineFastas:
    input: fasta=expand("{ref}/{name}", ref=REF, name=config["params"]["fasta"])
    output: "CombineFastas/database.fasta"
    run:
        fastas = input.fasta
        if not isinstance(input.fasta, list):
            fastas = [fastas]
        with open(str(output), 'w') as f:
            subprocess.check_call(['cat'] + fastas, stdout=f)


rule DecoyDatabase:
    input: "CombineFastas/database.fasta"
    output: "DecoyDatabase/database.fasta"
    params: params('DecoyDatabase')
    run:
        if config.get('fasta_has_decoy', False):
            shell("cp -- {} {}".format(input, output))
        else:
            openms.DecoyDatabase(input, output, ini=params)


rule MSGFPlusIndex:
    input: "DecoyDatabase/database.fasta"
    output: "DecoyDatabase/database.canno"
    shell:
        "java -Xmx3500M -cp ../usr/MSGFPlus.jar " +
            "edu.ucsd.msjava.msdbsearch.BuildSA -d {input}"


rule XTandemAdapter:
    input: fasta=rules.DecoyDatabase.output, mzml=data("{name}.mzML")
    output: "XTandemAdapter/{name}.idXML"
    params: params('XTandemAdapter')
    run:
        extra = ['-database', str(input.fasta)]
        key = "precursor_mass_tolerance"
        if key in params:
            extra += ["-" + key, str(params[key])]
        if 'xtandem_executable' in config:
            extra += ['-xtandem_executable', config['xtandem_executable']]
        openms.XTandemAdapter(input.mzml, output, extra_args=extra, ini=params)


rule MSGFPlusAdapter:
    input:
        fasta="DecoyDatabase/database.fasta",
        index="DecoyDatabase/database.canno",
        mzml=data("{name}.mzML")
    output: "MSGFPlusAdapter/{name}.idXML"
    params: params('MSGFPlusAdapter')
    run:
        extra = ['-database', str(input.fasta)]
        openms.MSGFPlusAdapter(
            input.mzml, output, extra_args=extra, ini=params
        )


rule IDPosteriorError:
    input: "{tool}/{name}.idXML"
    output: "IDPosteriorError_{tool}/{name}.idXML"
    params: params('IDPosteriorErrorProbability')
    run:
        openms.IDPosteriorErrorProbability(input, output, ini=params)


rule ConsensusID:
    input:
        tandem="IDPosteriorError_XTandemAdapter/{name}.idXML",
        msgf="IDPosteriorError_MSGFPlusAdapter/{name}.idXML"
    output: "ConsensusID/{name}.idXML"
    params: params('ConsensusID')
    run:
        openms.ConsensusID(input, output, ini=params)


rule IDFilter99:
    input: "ConsensusID/{name}.idXML"
    output: "IDFilter99/{name}.idXML"
    params: params('IDFilter99')
    run:
        openms.IDFilter(input, output, ini=params)


rule MapAlignerIdentification:
    input: expand("IDFilter99/{name}.idXML", name=INPUT_FILES)
    output:
        trafo="MapAligner/trafo.trafoXML",
        idXMLs=expand("MapAligner/{name}.idXML", name=INPUT_FILES)
    params: params('MapAlignerIdentification')
    run:
        extra = ['-trafo_out', str(output['trafo'])]
        openms.MapAlignerIdentification(
            input, output.idXMLs, extra_args=extra, ini=params
        )


rule MapRTTransformer:
    input:
        trafo="MapAligner/trafo.trafoXML",
        mzml=data("{name}.mzML")
    output: "MapRTTransformer/{name}.mzML"
    params: params('MapRTTransformer')
    run:
        extra = ['-trafo_in', input.trafo]
        openms.MapRTTransformer(
            input.mzml, output, extra_args=extra, ini=params
        )
        
        
rule MapRTTransformerID:
    input:
        trafo="MapAligner/trafo.trafoXML",
        idxml=expand("ConsensusID/{name}.idXML", name=INPUT_FILES)
    output: "MapRTTransformerID/{name}.idXML"
    params: params('MapRTTransformer')
    run:
        extra = ['-trafo_in', input.trafo]
        openms.MapRTTransformer(
            input.mzml, output, extra_args=extra, ini=params
        )


rule IDMerger_all:
    input: expand("MapRTTransformerID/{name}.idXML", name=INPUT_FILES)
    output: "IDMerger_all/all.idXML"
    params: params('IDMerger_all')
    run:
        openms.IDMerger(input, output, ini=params)


rule IDFilter:
    input: "IDMerger_all/all.idXML"
    output: "IDMerger_all/filtered.idXML"
    params: params('IDFilter60')
    run:
        openms.IDFilter(input, output, ini=params)


rule FeatureFinderIdentification:
    input:
        id="IDMerger_all/filtered.idXML",
        mzml="MapRTTransformer/{name}.mzML"
    output: "FeatureFinderID/{name}.featureXML"
    params: params('FeatureFinderIdentification')
    run:
        extra = ['-id', str(input.id)]
        openms.FeatureFinderIdentification(input.mzml, output, ini=params,
                                           extra_args=extra)


rule PeptideIndexer_all:
    input: fasta=rules.DecoyDatabase.output, \
           idxml="IDMerger_all/filtered.idXML"
    output: "PeptideIndexer/all.idXML"
    params: params('PeptideIndexer')
    run:
        extra=['-fasta', str(input.fasta)]
        openms.PeptideIndexer(input.idxml, output, extra_args=extra, ini=params)


rule Fido:
    input: "PeptideIndexer/all.idXML"
    output: result("proteins.idXML")
    params: params('FidoAdapter')
    run:
        openms.FidoAdapter(input, output, ini=params)


rule ProteinQuantifier:
    input:
        ms1=expand("FeatureFinderID/{name}.featureXML", name=INPUT_FILES),
        groups=result("proteins.idXML")
    output:
        proteins=result("proteins.csv"),
        peptides=result("peptides.csv")
    params: params('ProteinQuantifier')
    run:
        extra = [
            '-protein_groups', str(input['groups']),
            '-peptide_out', str(output['peptides'])
        ]
        openms.ProteinQuantifier(input, output, ini=params, extra_args=extra)


rule QCCalculator:
    input: mzml=os.path.join(DATA, "{name}.mzML"), \
           feature="IDMapper/{name}.featureXML", \
           idxml="IDFilter/{name}.idXML"
    output: "QCCalculator/{name}.qcML"
    params: params('QCCalculator')
    run:
        extra = ['-feature', str(input.feature), '-id', str(input.idxml)]
        openms.QCCalculator(input.mzml, output, extra_args=extra, ini=params)


def make_qc_plots(qcml, injection, run=None):
    accessions = {
        'features': 'QC:0000047',
        'identifications': 'QC:0000038',
        'mass_error': 'QC:0000038',
        'tic': 'QC:0000022',
        'raw_MS1': 'QC:0000044',
    }

    with tempfile.TemporaryDirectory() as tmp:
        for key, cv in accessions.items():
            csv = pjoin(tmp, key + '.csv')
            extra = ['-out_csv', csv, '-qp', cv]
            openms.QCExtractor(qcml, None, extra_args=extra)

        plot_dir = pjoin(tmp, 'plots')
        os.mkdir(plot_dir)
        shutil.copy(injection, pjoin(tmp, 'injection_times.csv'))

        for file in os.listdir(R_HOME):
            if os.path.splitext(file)[1].lower() == '.r':
                script = os.path.join(R_HOME, file)
                try:
                    subprocess.check_output(['Rscript', script, tmp, plot_dir])
                except subprocess.CalledProcessError:
                    print("Rscript returned non-zero. Ignoring.")

        plots = {}
        files = os.listdir(plot_dir)
        for name in sorted(files):
            stem, ext = os.path.splitext(name)
            path = os.path.join(plot_dir, name)
            plot = plots.setdefault(stem, {})
            with open(path, 'rb') as f:
                if ext.lower() == '.png':
                    png = f.read()
                    plot['png'] = base64.encodebytes(png).decode()
                elif ext.lower() == '.svg':
                    plot['svg'] = f.read().decode()
                elif ext.lower() in ['.txt', '.html']:
                    plot['desc'] = f.read().decode()
                elif ext.lower() == '.csv':
                    plot['table'] = list(csv.reader(f))
                elif ext.lower() == '.title':
                    plot['title'] = f.read().decode()

    if any("title" not in val for val in plots.values()):
        raise ValueError("Invalid output of R script. Title is missing")
    return plots


rule FixQCML:
    input: "QCCalculator/{name}.qcML"
    output: "QCCalculator_fixed/{name}.qcML"
    shell:
        'grep -Fv "UTF-8" {input} > {output}'


rule HTML:
    input: qcml="QCCalculator_fixed/{name}.qcML", \
           injection="InjectionTime/{name}.csv"
    output: os.path.join(RESULT, "{name}.html")
    run:
        import jinja2
        NAMESPACE = "{http://www.prime-xs.eu/ms/qcml}"  # openms 1.12?
        #NAMESPACE = ""
        tree = ElementTree.ElementTree(file=str(input['qcml']))
        runs = tree.findall(NAMESPACE + 'runQuality')

        fastas = expand("{ref}/{name}", name=config['params']['fasta'], ref=REF)
        if not isinstance(fastas, list):
            fastas = list(fastas)

        fasta_md5 = []
        fasta_size = []
        for fasta in fastas:
            md5 = subprocess.check_output(['md5sum', fasta])
            fasta_md5.append(md5.split()[0].decode())
            size = os.stat(fasta).st_size / 1024 / 1024
            fasta_size.append("%.2fM" % size)

        qcprot = {
            'date': datetime.strftime(datetime.now(), "%d. %B %Y at %H:%M"),
            'version': VERSION,
            'fasta_names': fastas,
            'fasta_md5s': fasta_md5,
            'fasta_sizes': fasta_size,
        }

        orig_ini = pjoin(R_HOME, '..', 'inis')
        popen = subprocess.Popen(
            ['diff', '-u', '-w', orig_ini, INI_PATH],
            stdout=subprocess.PIPE
        )
        ini_diff = popen.communicate()[0].decode()
        if ini_diff.strip() != "":

            p = subprocess.Popen(
                ['pygmentize', '-l', 'diff', '-f', 'html'],
                stdin=subprocess.PIPE, stdout=subprocess.PIPE
            )
            out, _ = p.communicate(ini_diff.encode())
            if p.returncode == 0:
                qcprot['ini_diff'] = out.decode()
            else:
                qcprot['ini_diff'] = ini_diff

        p = subprocess.Popen(
            ['pygmentize', '-l', 'json', '-f', 'html'],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE
        )
        out, _ = p.communicate(json.dumps(config, indent=4).encode())
        qcprot['config'] = out.decode()

        qcprot['runs'] = []

        for run_el in runs:  # openms 1.12?:
            #for run_el in tree.findall('RunQuality'):
            run = {
                'id': run_el.get('ID'),
                'quality_params': [],
            }
            qcprot['runs'].append(run)
            run['plots'] = make_qc_plots(str(input['qcml']),
                                         str(input['injection']),
                                         run=run['id'])
            data = run['quality_params']

            #openms 1.12
            quality_params = run_el.findall(NAMESPACE + 'qualityParameter')
            #quality_params = run_el.findall('QualityParameter')
            for param_el in quality_params:
                if 'value' in param_el.keys() and 'name' in param_el.keys():
                    data.append((param_el.get('name'), param_el.get('value')))

        with open(pjoin(R_HOME, 'report_template.html')) as f:
            template = jinja2.Template(f.read())
            with open(str(output), 'w') as f_out:
                f_out.write(template.render(qcprot=qcprot))
