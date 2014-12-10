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

configfile: "config.json"

for key in ['R_HOME', 'OPENMS_BIN', 'QCPROT_VERSION']:
    if key not in os.environ:
        print("Environment variable %s is not defined. Aborting." % key,
              file=sys.stderr)
        exit(1)

if not pexists('work'):
    os.mkdir("work")

R_HOME = os.environ['R_HOME']
OPENMS_BIN = os.environ['OPENMS_BIN']
QCPROT_VERSION = os.environ['QCPROT_VERSION']

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

INI_PATH = config.get('ini_path', os.environ.get('INI_PATH', 'inis'))
openms = OpenMS(OPENMS_BIN, INI_PATH, 'logs')

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
for name in os.listdir('mzml'):
    if name.lower().endswith('.mzml'):
        INPUT_FILES.append(os.path.basename(name)[:-5])
        if not name.endswith('.mzML'):
            print("Extension mzML is case sensitive.", file=sys.stderr)
            exit(1)


rule all:
    input: ["result/{name}.html".format(name=name) for name in INPUT_FILES]


rule FileFilter:
    input: "mzml/{name}.mzML"
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


rule IDPosteriorError:
    input: "work/XTandemAdapter/{name}.idXML"
    output: "work/IDPosteriorXTandem/{name}.idXML"
    params: params('IDPosteriorErrorProbability')
    run:
        openms.IDPosteriorErrorProbability(input, output, ini=params)


rule PeptideIndexer:
    input: fasta=rules.DecoyDatabase.output, \
           idxml="work/IDPosteriorXTandem/{name}.idXML"
    output: "work/PeptideIndexer/{name}.idXML"
    params: params('PeptideIndexer')
    run:
        extra=['-fasta', str(input.fasta)]
        openms.PeptideIndexer(input.idxml, output, extra_args=extra, ini=params)


rule FalseDiscoveryRate:
    input: "work/PeptideIndexer/{name}.idXML"
    output: "work/FalseDiscoveryRate/{name}.idXML"
    params: params('FalseDiscoveryRate')
    run:
        openms.FalseDiscoveryRate(input, output, ini=params)


rule IDFilter:
    input: "work/FalseDiscoveryRate/{name}.idXML"
    output: "work/IDFilter/{name}.idXML"
    params: params('IDFilter')
    run:
        openms.IDFilter(input, output, ini=params)


rule IDMapper:
    input: feature="work/FeatureFinderCentroided/{name}.featureXML", \
           idxml="work/IDFilter/{name}.idXML"
    output: 'work/IDMapper/{name}.featureXML'
    params: params('IDMapper')
    run:
        extra = ['-id', str(input.idxml)]
        openms.IDMapper(input.feature, output, extra_args=extra, ini=params)


rule QCCalculator:
    input: mzml="mzml/{name}.mzML", \
           feature="work/IDMapper/{name}.featureXML", \
           idxml="work/IDFilter/{name}.idXML"
    output: "work/QCCalculator/{name}.qcML"
    params: params('QCCalculator')
    run:
        extra = ['-feature', input.feature, '-id', input.idxml]
        openms.QCCalculator(input.mzml, output, extra_args=extra, ini=params)


def make_qc_plots(qcml, run=None):
    accessions = {
        'fractional_mass': 'QC:0000047',
        'mass_acc': 'QC:0000038',
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
                #subprocess.check_call(['Rscript', file, tmp, plot_dir])
                pass

        plots = {}
        files = os.listdir(plot_dir)
        for file in sorted(files):
            path, name = os.path.split(file)
            stem, ext = os.path.splitext(name)
            plot = plots.setdefault(stem, {})
            with open(file) as f:
                if ext.lower() == '.png':
                    png = f.read().encode()
                    plot['png'] = base64.encodbytes(png).decode()
                elif ext.lower() == '.svg':
                    plot['svg'] = f.read()
                elif ext.lower() in ['.txt', '.html']:
                    plot['desc'] = f.read()
                elif ext.lower() == '.csv':
                    plot['table'] = list(csv.reader(f))
                elif ext.lower() == '.title':
                    plot['title'] = f.read()

    return plots


rule HTML:
    input: "work/QCCalculator/{name}.qcML"
    output: "result/{name}.html"
    run:
        NAMESPACE = "{http://www.prime-xs.eu/ms/qcml}"
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

        qcprot['runs'] = []

        for run_el in runs:
            run = {
                'id': run_el.get('ID'),
                'quality_params': [],
            }
            qcprot['runs'].append(run)
            run['plots'] = make_qc_plots(str(input), run=run['id'])
            data = run['quality_params']
            quality_params = run_el.findall(NAMESPACE + 'qualityParameter')
            for param_el in quality_params:
                if 'value' in param_el.keys() and 'name' in param_el.keys():
                    data.append((param_el.get('name'), param_el.get('value')))

        with open(pjoin(R_HOME, 'report_template.html')) as f:
            template = jinja2.Template(f.read())
            with open(str(output), 'w') as f_out:
                print(qcprot)
                f_out.write(template.render(qcprot=qcprot))
