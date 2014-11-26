import json
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

configfile: "config.json"

for key in ['R_HOME', 'OPENMS_BIN', 'QCPROT_VERSION']:
    if key not in os.environ:
        print("Environment variable %s is not defined. Aborting." % key,
              file=sys.stderr)
        exit(1)

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
                command += ['-ini', json.loads(str(ini))['path']]
            command += extra_args

            log_std = pjoin(self._log_dir, logfile + '.out')
            log_err = pjoin(self._log_dir, logfile + '.err')
            with open(log_std, 'w') as out:
                with open(log_err, 'w') as err:
                    subprocess.check_call(command, stdout=out, stderr=err)

        return wrapper

INI_PATH = config.get('ini_path', os.environ('INI_PATH', 'inis'))
openms = OpenMS(OPENMS_BIN, INI_PATH, 'logs')

# store the content of the ini file, so that snakemake will run
# rules agrain if the parameters inside the file change.
def params(name):
    path = pjoin(INI_PATH, name + '.ini')
    try:
        with open(path, 'r') as f:
            # TODO replace makes sure that there are no wildcards
            # in the file content. This is not exactly clean ;-)
            data = {'path': path, 'content': f.read().replace('{', '[')}
            return json.dumps(data, sort_keys=True)
    except FileNotFoundError as e:
        raise ValueError("ini file '%s' not found" % path) from e


INPUT_FILES = []
for name in os.listdir('mzml'):
    if name.lower().endswith('.mzml'):
        INPUT_FILES.append(os.path.basename(name)[:-5])


rule all:
    input: ["result/{name}.html".format(name=name) for name in INPUT_FILES]


rule FeatureFinderCentroided:
    input: "mzml/{name}.mzml"
    output: "FeatureFinderCentroided/{name}.featureXML"
    params: params('FeatureFinderCentroided')
    run:
        openms.FeatureFinderCentroided(input, output, ini=params)


rule DecoyDatabase:
    input: config["fasta"]
    output: "DecoyDatabase/database.fasta"
    params: params('DecoyDatabase')
    run:
        if config.get('fasta_has_decoy', False):
            shell("cp -- {} {}".format(input, output))
        else:
            openms.DecoyDatabase(input, output, ini=params)


rule XTandemAdapter:
    input: fasta=rules.DecoyDatabase.output, mzml="mzml/{name}.mzml"
    output: "XTandemAdapter/{name}.idXML"
    params: params('XTandemAdapter')
    run:
        extra = ['-database', str(input.fasta)]
        openms.XTandemAdapter(input.mzml, output, extra_args=extra, ini=params)


rule IDPosteriorError:
    input: "XTandemAdapter/{name}.idXML"
    output: "IDPosteriorXTandem/{name}.idXML"
    params: params('IDPosteriorErrorProbability')
    run:
        openms.IDPosteriorErrorProbability(input, output, ini=params)


rule PeptideIndexer:
    input: fasta=rules.DecoyDatabase.output, \
           idxml="IDPosteriorXTandem/{name}.idXML"
    output: "PeptideIndexer/{name}.idXML"
    params: params('PeptideIndexer')
    run:
        extra=['-fasta', str(input.fasta)]
        openms.PeptideIndexer(input.idxml, output, extra_args=extra, ini=params)


rule FalseDiscoveryRate:
    input: "PeptideIndexer/{name}.idXML"
    output: "FalseDiscoveryRate/{name}.idXML"
    params: params('PeptideIndexer')
    run:
        openms.FalseDiscoveryRate(input, output, ini=params)


rule IDFilter:
    input: "FalseDiscoveryRate/{name}.idXML"
    output: "IDFilter/{name}.idXML"
    params: params('IDFilter')
    run:
        openms.IDFilter(input, output, ini=params)


rule IDMapper:
    input: feature="FeatureFinderCentroided/{name}.featureXML", \
           idxml="IDFilter/{name}.idXML"
    output: 'IDMapper/{name}.featureXML'
    params: params('IDMapper')
    run:
        extra = ['-id', str(input.idxml)]
        openms.IDMapper(input.feature, output, extra_args=extra, ini=params)


rule QCCalculator:
    input: mzml="mzml/{name}.mzml", \
           feature="IDMapper/{name}.featureXML", \
           idxml="IDFilter/{name}.idXML"
    output: "QCCalculator/{name}.qcML"
    params: params('QCCalculator')
    run:
        extra = ['-feature', input.feature, '-id', input.idxml]
        openms.QCCalculator(input.mzml, output, extra_args=extra, ini=params)


rule PlotQC:
    input: "QCCalculator/{name}.qcML"
    output: "result/{name}.qcML"
    run:
        extra_cv = {
            'fractional_mass':
                ('QC:0000047', 'fractionalMass.R'),
            'mass_acc':
                ('QC:0000044', 'MassAcc.R'),
            'mass_error':
                ('QC:0000038', 'MassError.R'),
            'tic':
                ('QC:0000022', 'TIC.R')
        }

        with tempfile.TemporaryDirectory() as tmp:
            for key, (cv, script) in extra_cv.items():
                csv = pjoin(tmp, key + '.csv')
                png = pjoin(tmp, key + '.png')
                extra = ['-out_csv', csv, '-qp', cv]
                openms.QCExtractor(input, None, extra_args=extra)
                script_path = pjoin(R_HOME, script)
                subprocess.check_call(['Rscript', script_path, csv, png])

            script_path = pjoin(R_HOME, 'IDRatio.R')
            command = [
                'Rscript',
                script_path,
                pjoin(tmp, 'mass_acc.csv'),
                pjoin(tmp, 'mass_error.csv'),
                pjoin(tmp, 'id_ratio.png'),
            ]
            subprocess.check_call(command)

            tmp1_qcml = pjoin(tmp, 'tmp1.qcml')
            tmp2_qcml = pjoin(tmp, 'tmp2_qcml')
            shutil.copy(str(input), tmp1_qcml)
            args = [
                ('tic.png', 'MS:1000235', 'QC:0000023'),
                ('mass_error.png', 'QC:0000054', 'QC:0000041'),
                ('mass_acc.png', 'QC:0000053', 'QC:0000041'),
                ('fractional_mass.png', 'QC:0000043', 'QC:0000007'),
                ('id_ratio.png', 'QC:0000053', 'QC:0000035')
            ]
            for png, cv_acc, qp_att_acc in args:
                png_path = pjoin(tmp, png)
                extra = ['-plot', png_path,
                         '-cv_acc', cv_acc,
                         '-qp_att_acc', qp_att_acc]
                openms.QCEmbedder(tmp1_qcml, tmp2_qcml, extra_args=extra)
                shutil.copy(tmp2_qcml, tmp1_qcml)

            openms.QCShrinker(tmp1_qcml, output)


rule HTML:
    input: "result/{name}.qcML"
    output: "result/{name}.html"
    run:
        tree = ElementTree(file=str(input))
        runs = tree.findall('RunQuality')

        fasta_md5 = subprocess.check_output(['md5sum', config['fasta']])
        fasta_size = os.stat(config['fasta']).st_size
        qcprot = {
            'date': datetime.strftime(datetime.now(), "%d. %B %Y at %H:%M"),
            'version': QCPROT_VERSION,
            'fasta_name': config['fasta'],
            'fasta_md5': fasta_md5.split()[0],
            'fasta_size': "%.2fM" % (fasta_size / 1024 / 1024),
        }
        qcprot['runs'] = []
        for run_el in runs:
            run = {
                'id': run_el.get('ID'),
                'data': [],
                'attachments': [],
            }
            qcprot['runs'].append(run)
            data = run['data']
            quality_params = run_el.findall('QualityParameter')
            attachments = run_el.findall('Attachment')
            for param_el in quality_params:
                if 'value' in param_el.keys() and 'name' in param_el.keys():
                    data.append((param_el.get('name'), param_el.get('value')))
            for attach_el in attachments:
                attach = {}
                run['attachments'].append(attach)
                attach['name'] = attach_el.get('name')
                binary = attach_el.find('binary')
                if binary is not None:
                    attach['plot_data'] = binary.text.strip()
                columns = attach_el.find('TableColumnTypes')
                if columns is not None:
                    table = []
                    attach['table'] = table
                    table.append(columns.text.split())
                    for row in attach_el.findall('TableRowValues'):
                        table.append(row.text.split())

        with open('report_template.html') as f:
            template = jinja2.Template(f.read())
            with open(str(output), 'w') as f_out:
                f_out.write(template.render(qcprot=qcprot))
