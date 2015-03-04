import os
import sys
import subprocess
import tempfile
import uuid
import shutil
import jinja2
from datetime import datetime
from os.path import join as pjoin
from os.path import exists as pexists
from xml.etree.ElementTree import ElementTree
import hashlib
import base64
import csv
import glob
import functools

BASE_DIR = os.getcwd()

configfile: "config.json"
workdir: "../../var/qcprot"


if "OPENMS_BIN" not in os.environ:
    path = subprocess.check_output(["which", "IDMerger"]).decode()
    OPENMS_BIN = os.path.split(path)[0]
else:
    OPENMS_BIN = os.environ['OPENMS_BIN']

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
            if output is not None and not isinstance(output, list):
                if not pexists(os.path.dirname(str(output))):
                    os.mkdir(os.path.dirname(str(output)))
            elif output is not None:
                for out in output:
                    if not pexists(os.path.dirname(out)):
                        os.mkdir(os.path.dirname(out))
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
                if isinstance(input, list):
                    command += ['-in'] + list(input)
                else:
                    command += ['-in'] + [str(input)]
            if output is not None:
                if isinstance(output, list):
                    command += ['-out'] + list(output)
                else:
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
                    print(command)
                    subprocess.check_call(command, stdout=out, stderr=err)

        return wrapper

INI_PATH = os.path.join(BASE_DIR, config.get('ini_path', os.environ.get('INI_PATH', 'inis')))
RESULTS = "../../results/qcprot"

if not os.path.exists(RESULTS):
    os.mkdir(RESULTS)

openms = OpenMS(OPENMS_BIN, INI_PATH, os.path.join(RESULTS, 'logs'))

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


DATA = '../../data'


INPUT_FILES = []
for name in os.listdir(DATA):
    if name.lower().endswith('.mzml'):
        INPUT_FILES.append(os.path.basename(name)[:-5])
        if not name.endswith('.mzML'):
            print("Extension mzML is case sensitive.", file=sys.stderr)
            exit(1)


samples = config['samples']


rule all:
    input: expand(os.path.join(RESULTS, "{sample}.csv"), sample=samples.keys()), \
           "work/FalseDiscoveryRateAll/all.idXML"


rule IdToResults:
    input: "work/{target}/{sample}.idXML"
    output: os.path.join(RESULTS, "{sample}_{target}.idXML")
    shell: "ln {input} {output}"


rule FileFilter:
    input: os.path.join(DATA, "{name}.mzML")
    output: "work/FileFilter/{name}.mzML"
    run:
        openms.FileFilter(input, output, extra_args=['-sort'])


rule PeakPicker:
    input: "work/FileFilter/{name}.mzML"
    output: "work/PeakPicker/{name}.mzML"
    params: params('PeakPickerHiRes')
    run:
        openms.PeakPickerHiRes(input, output, ini=params)


rule FeatureFinderCentroided:
    input: "work/PeakPicker/{name}.mzML"
    output: "work/FeatureFinderCentroided/{name}.featureXML"
    params: params('FeatureFinderCentroided')
    run:
        openms.FeatureFinderCentroided(input, output, ini=params)


rule CombineFastas:
    input: fasta=config["fasta"]
    output: "work/CombineFastas/database.fasta"
    run:
        fastas = input.fasta
        if not isinstance(input.fasta, list):
            fastas = [fastas]
        with open(str(output), 'w') as f:
            subprocess.check_call(['cat'] + fastas, stdout=f)


rule DecoyDatabase:
    input: "work/CombineFastas/database.fasta"
    output: "work/DecoyDatabase/database.fasta"
    params: params('DecoyDatabase')
    run:
        if config.get('fasta_has_decoy', False):
            shell("cp -- {} {}".format(input, output))
        else:
            openms.DecoyDatabase(input, output, ini=params)


rule XTandemAdapter:
    input: fasta=rules.DecoyDatabase.output, mzml="work/PeakPicker/{name}.mzML"
    output: "work/XTandemAdapter/{name}.idXML"
    params: params('XTandemAdapter')
    run:
        extra = ['-database', str(input.fasta)]
        if 'xtandem_executable' in config:
            extra += ['-xtandem_executable', config['xtandem_executable']]
        openms.XTandemAdapter(input.mzml, output, extra_args=extra, ini=params)


def files_by_sample(format_string):
    def inner(wildcards):
        if isinstance(wildcards, str):
            sample = wildcards
        else:
            sample = wildcards['sample']
        names = []
        for pattern in config["samples"][sample]:
            for file in glob.glob(os.path.join(DATA, pattern)):
                path, name = os.path.split(file)
                base, ext = os.path.splitext(name)
                names.append(base)
        r = expand(format_string, name=names)
        return r
    return inner


def all_files(format_string):
    def inner(wildcards):
        samples = config['samples'].keys()
        print(samples)
        lists = (files_by_sample(format_string)(sample) for sample in samples)
        r = sum(lists, [])
        print(r)
        return r
    return inner


rule IDMerger:
    input: files_by_sample("work/XTandemAdapter/{name}.idXML")
    output: "work/IDMerger/{sample}.idXML"
    params: params("IDMerger")
    run:
        openms.IDMerger(list(input), output, ini=params)


rule IDMergerAll:
    input: all_files("work/XTandemAdapter/{name}.idXML")
    output: "work/IDMergerAll/all.idXML"
    params: params("IDMerger")
    run:
        openms.IDMerger(list(input), output, ini=params)


rule PeptideIndexerAll:
    input: fasta=rules.DecoyDatabase.output, \
           idxml="work/IDMergerAll/all.idXML"
    output: "work/PeptideIndexerAll/all.idXML"
    params: params('PeptideIndexer')
    run:
        extra=['-fasta', str(input.fasta)]
        openms.PeptideIndexer(input.idxml, output, extra_args=extra, ini=params)


rule FalseDiscoveryRateAll:
    input: "work/PeptideIndexerAll/all.idXML"
    output: "work/FalseDiscoveryRateAll/all.idXML"
    params: params('FalseDiscoveryRate')
    run:
        openms.FalseDiscoveryRate(input, output, ini=params)


rule PeptideIndexer:
    input: fasta=rules.DecoyDatabase.output, \
           idxml="work/IDMerger/{sample}.idXML"
    output: "work/PeptideIndexer/{sample}.idXML"
    params: params('PeptideIndexer')
    run:
        extra=['-fasta', str(input.fasta)]
        openms.PeptideIndexer(input.idxml, output, extra_args=extra, ini=params)


rule FalseDiscoveryRate:
    input: "work/PeptideIndexer/{sample}.idXML"
    output: "work/FalseDiscoveryRate/{sample}.idXML"
    params: params('FalseDiscoveryRate')
    run:
        openms.FalseDiscoveryRate(input, output, ini=params)


rule IDRipper:
    input: "work/FalseDiscoveryRate/{sample}.idXML"
    output: "work/IDRipper/{sample}"
    params: params('IDRipper')
    run:
        os.mkdir(output[0])
        extra = ['-out_path', output[0] + '/bug_ignored']
        openms.IDRipper(input, None, ini=params, extra_args=extra)
        outfiles = files_by_sample(output[0] + '/{name}.idXML')(wildcards)
        for file in outfiles:
            assert os.path.exists(file), "Not found: %s" % file


rule IDMapper:
    input: "work/IDRipper/{sample}", \
           files_by_sample("work/FeatureFinderCentroided/{name}.featureXML")
    output: 'work/IDMapper/{sample}'
    params: params('IDMapper')
    run:
        idfiles = files_by_sample(input[0] + "/{name}.idXML")(wildcards)
        features = files_by_sample(
                "work/FeatureFinderCentroided/{name}.featureXML")(wildcards)
        outfiles = files_by_sample(output[0] + "/{name}.featureXML")(wildcards)
        for file in idfiles + features:
            assert os.path.exists(file), "File not found: %s" % file
        for idfile, feature, out in zip(idfiles, features, outfiles):
            extra = ['-id', str(idfile)]
            openms.IDMapper(feature, out, extra_args=extra, ini=params)
        for outfile in outfiles:
            assert os.path.exists(outfile), \
                    "IDMapper did not produce file %s" % outfile


rule MapAligner:
    input: "work/IDMapper/{sample}"
    output: "work/MapAligner/{sample}"
    params: params('MapAlignerPoseClustering')
    run:
        infiles = files_by_sample(input[0] + "/{name}.featureXML")(wildcards)
        outfiles = files_by_sample(output[0] + "/{name}.featureXML")(wildcards)
        for file in infiles:
            assert os.path.exists(file), "File not found: %s" % file

        openms.MapAlignerPoseClustering(infiles, outfiles, ini=params)
        for file in outfiles:
            assert os.path.exists(file), "File not found: %s" % file


rule FeatureLinker:
    input: "work/MapAligner/{sample}"
    output: "work/FeatureLinker/{sample}.consensusXML"
    params: params('FeatureLinkerUnlabledQT')
    run:
        infiles = files_by_sample(input[0] + "/{name}.featureXML")(wildcards)
        for file in infiles:
            assert os.path.exists(file), "File not found: %s" % file

        openms.FeatureLinkerUnlabeledQT(infiles, output, ini=params)


rule ConsensusMapNormalizer:
    input: "work/FeatureLinker/{sample}.consensusXML"
    output: "work/ConsensusMapNormalizer/{sample}.consensusXML"
    params: params('ConsensusMapNormalizer')
    run:
        openms.ConsensusMapNormalizer(input, output, ini=params)


rule TextExporter:
    input: "work/ConsensusMapNormalizer/{sample}.consensusXML"
    output: os.path.join(RESULTS, "{sample}.csv")
    params: params('TextExporter')
    run:
        openms.TextExporter(input, output, ini=params)


rule QCCalculator:
    input: mzml=os.path.join(DATA, "{name}.mzML"), \
           feature="work/IDMapper/{name}.featureXML", \
           idxml="work/IDFilter/{name}.idXML"
    output: "work/QCCalculator/{name}.qcML"
    params: params('QCCalculator')
    run:
        extra = ['-feature', input.feature, '-id', input.idxml]
        openms.QCCalculator(input.mzml, output, extra_args=extra, ini=params)


def make_qc_plots(qcml, run=None):
    accessions = {
        'features': 'QC:0000047',
        'identifications': 'QC:0000038',
        'mass_error': 'QC:0000038',
        'tic': 'QC:0000022',
    }

    with tempfile.TemporaryDirectory() as tmp:
        for key, cv in accessions.items():
            csv = pjoin(tmp, key + '.csv')
            extra = ['-out_csv', csv, '-qp', cv]
            openms.QCExtractor(qcml, None, extra_args=extra)

        plot_dir = pjoin(tmp, 'plots')
        os.mkdir(plot_dir)

        for file in os.listdir(R_HOME):
            if os.path.splitext(file)[1].lower() == '.r':
                script = os.path.join(R_HOME, file)
                subprocess.call(['Rscript', script, tmp, plot_dir])

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


rule HTML:
    input: "work/QCCalculator/{name}.qcML"
    output: os.path.join(RESULTS, "{name}.html")
    run:
        #NAMESPACE = "{http://www.prime-xs.eu/ms/qcml}"  # openms 1.12?
        NAMESPACE = ""
        tree = ElementTree(file=str(input))
        runs = tree.findall(NAMESPACE + 'runQuality')

        fastas = config['fasta']
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
            'version': QCPROT_VERSION,
            'fasta_names': fastas,
            'fasta_md5s': fasta_md5,
            'fasta_sizes': fasta_size,
        }

        orig_ini = pjoin(R_HOME, '..', 'inis')
        ini_diff = subprocess.check_output(
            ['diff', '-u', '-w', orig_ini, INI_PATH]
        ).decode()
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

        #for run_el in runs:  # openms 1.12?
        for run_el in tree.findall('RunQuality'):
            run = {
                'id': run_el.get('ID'),
                'quality_params': [],
            }
            qcprot['runs'].append(run)
            run['plots'] = make_qc_plots(str(input), run=run['id'])
            data = run['quality_params']

            #openms 1.12
            #quality_params = run_el.findall(NAMESPACE + 'qualityParameter')
            quality_params = run_el.findall('QualityParameter')
            for param_el in quality_params:
                if 'value' in param_el.keys() and 'name' in param_el.keys():
                    data.append((param_el.get('name'), param_el.get('value')))

        with open(pjoin(R_HOME, 'report_template.html')) as f:
            template = jinja2.Template(f.read())
            with open(str(output), 'w') as f_out:
                f_out.write(template.render(qcprot=qcprot))
