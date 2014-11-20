import json
import os


with open("config.json") as f:
    CONFIG = json.load(f)

config['R_home'] = '/lustre_cfc/software/qbic/qcprot/r_scripts'

def file_name_and_ext(path):
    head, tail = os.path.split(path)
    if not tail:
        return (None, None) 
    return os.path.splitext(tail)

def file_name_w_ext(path):
    head, tail = os.path.split(path)
    if not tail:
        return None
    return os.path.splitext(tail)[0]

def createOutputFileNames(path, input_files, extension):
    lis = []
    for infi in input_files:
        lis.append(path + file_name_w_ext(infi)+ extension)
    return lis

"""This function merges the output of different search engines for the same initial file. It assumes, that they were written in the same order"""
def mergeIDPepOutput(lis):
    if (not isinstance(lis, list)) or lis == []: 
         return None
    #the idea was that all files have the same name. However, that does not need to be true i guess?
    #for l in lis:
    #    sort(l)
    ret_lis = []
    length = len(lis[1])
    for i in range(length):
        new_line = []
        for j in range(len(lis)):
            new_line.append(lis[j][i])
        ret_lis.append(new_line)
    return ret_lis


INPUT_FILES = CONFIG["input_files"]
RESULT_FILES = createOutputFileNames("results/",INPUT_FILES,".qcML")
RESULT_FILE = "results/qcmlresult.idXML"
#prepended module load to all shell commands 
shell.prefix(CONFIG["shell_prefix"] + "; " )

rule all:
    input: RESULT_FILES

FFCOutput = createOutputFileNames("FFC/",INPUT_FILES, ".featureXML")
rule FeatureFinderCentroided:
    input: INPUT_FILES
    output: FFCOutput
    run:
        for inp, outp in zip (input, output):
            shell("FeatureFinderCentroided -in {inp} -out {outp}")

rule DecoyDatabase:
    input: CONFIG["REF"]
    output: "DecoyDatabase/%s%s" % file_name_and_ext(CONFIG["REF"])#"DecoyDatabase/%s_decoy%s" % file_name_and_ext(CONFIG["REF"])
    shell: "cp {input}* DecoyDatabase/"#"DecoyDatabase -in {input} -out {output}"

XTandemOutput = createOutputFileNames("XTandemAdapter/",INPUT_FILES,".idXML")

rule XTandemAdapter:
    input: INPUT_FILES, REF = rules.DecoyDatabase.output
    output: XTandemOutput
    threads: 1
    run:
        for inp, outp in zip(input,output):
            shell("{rules.XTandemAdapter.name} -in {inp} -out {outp} -database {input.REF} -threads {threads}")

"""
MyriMatchOutput = createOutputFileNames("MyriMatchAdapter/", INPUT_FILES, ".idXML")

rule MyrimatchAdapter:
    input: INPUT_FILES, REF = rules.DecoyDatabase.output
    output: MyriMatchOutput
    threads: 1
    run:
        for inp, outp in zip(input,output):
            shell("MyriMatchAdapter -in {inp} -out {outp} -database {input.REF}")
"""

IDPepOfXtandemOutput = createOutputFileNames("IDPepOfXTandem/", XTandemOutput,".idXML")

rule IDPEPOfXTandem:
    input: XTandemOutput
    output: IDPepOfXtandemOutput
    run:
        for inp, outp in zip(input,output):
            shell("IDPosteriorErrorProbability -in {inp} -out {outp} -output_name IDPepOfXTandem/delete.txt -threads {threads}")

"""
IDPepOfMyriMatchOutput = createOutputFileNames("IDPepOfMyriMatch/", MyriMatchOutput,".idXML")

rule IDPEPOfMyriMatch:
    input: MyriMatchOutput
    output: IDPepOfMyriMatchOutput
    run:
        for inp, outp in zip(input,output):
            shell("IDPosteriorErrorProbability -in {inp} -out {outp} -output_name IDPepOfMyriMatch/delete.txt -threads {threads}")

IDMergerOutput = createOutputFileNames("IDMerger/",INPUT_FILES, ".idXML")
rule IDMerger:
    input: mergeIDPepOutput([IDPepOfXtandemOutput, IDPepOfMyriMatchOuput])
    output: IDMergerOutput
    threads: 1
    run:
        for inp,outp in zip(input,output):
            shell("IDMerger -in {inp} -out {outp} -threads {threads}")
"""

ConsensusIDOutput = createOutputFileNames("ConsensusID/",IDPepOfXtandemOutput, ".idXML")
rule ConsensusID:
    input: IDPepOfXtandemOutput
    output: ConsensusIDOutput
    threads: 1
    run:
        for inp, outp in zip(input, output):
            shell("ConsensusID -in {inp} -out {outp}")

PepIndexerOutput = createOutputFileNames("PepIndexer/",ConsensusIDOutput,".idXMl")
rule PeptideIndexer:
    input: ConsensusIDOutput
    output: PepIndexerOutput
    run:
        for inp, outp in zip(input, output):
            shell("PeptideIndexer -in {inp} -out {outp} -fasta {rules.DecoyDatabase.output} -decoy_string sw -prefix -enzyme:specificity none -missing_decoy_action warn")

FDROutput = createOutputFileNames("FalseDiscoveryRate/",PepIndexerOutput,".idXML")
rule FDR:
    input: PepIndexerOutput
    output: FDROutput
    run:
        for inp, outp in zip(input, output):
            shell("FalseDiscoveryRate -in {inp} -out {outp}")

IDFilterOutput = createOutputFileNames("IDFilter/",PepIndexerOutput, ".idXML")
rule IDFilter:
    input: FDROutput 
    output: IDFilterOutput
    run:
        for inp, outp in zip(input, output):
            shell("IDFilter -in {inp} -out {outp}")

IDMapperOutput = createOutputFileNames("IDMapper/",FFCOutput, ".featureXML")
rule IDMapper:
    input: _in = FFCOutput, _id=IDFilterOutput
    output: IDMapperOutput
    run:
        for _id,_in, outp in zip(input._id, input._in,output):
            shell("IDMapper -id {_id} -in {_in} -out {outp}")


QCCalculatorOutput = createOutputFileNames("QCCalculator/",IDMapperOutput, ".qcML")
rule QCCalculator:
    input: _in = INPUT_FILES, _feature = IDMapperOutput, _id = IDFilterOutput
    output: QCCalculatorOutput
    run:
        for _in, _feature,_id, outp in zip(input._in, input._feature,input._id, output):
            shell("QCCalculator -in {_in} -feature {_feature} -id {_id} -out {outp}")

QCExtractorFractionalMassOutput = createOutputFileNames("QCExtractorFractionalMass/", QCCalculatorOutput, ".csv")
rule QCExtractorFractionalMass:
    input: QCCalculatorOutput
    output: QCExtractorFractionalMassOutput
    run:
        for inp, outp in zip(input, output):
            shell("QCExtractor -in {inp} -out_csv {outp} -qp QC:0000047")

#Mass Acc
QCExtractorIDR44Output = createOutputFileNames("QCExtractorIDR44/", QCCalculatorOutput, ".csv")
rule QCExtractorIDR44:
    input: QCCalculatorOutput
    output: QCExtractorIDR44Output
    run:
        for inp, outp in zip(input, output):
            shell("QCExtractor -in {inp} -out_csv {outp} -qp QC:0000044")

#Mass Error
QCExtractorIDR38Output = createOutputFileNames("QCExtractorIDR38/", QCCalculatorOutput, ".csv")
rule QCExtractorIDR38:
    input: QCCalculatorOutput
    output: QCExtractorIDR38Output
    run:
        for inp, outp in zip(input, output):
            shell("QCExtractor -in {inp} -out_csv {outp} -qp QC:0000038")

QCExtractorTICOutput = createOutputFileNames("QCExtractorTIC/", QCCalculatorOutput, ".csv")
rule QCExtractorTIC:
    input: QCCalculatorOutput
    output: QCExtractorTICOutput
    run:
        for inp, outp in zip(input, output):
            shell("QCExtractor -in {inp} -out_csv {outp} -qp QC:0000022")


fractionalMassOutput = createOutputFileNames("fractionalMass/",QCExtractorFractionalMassOutput,".png")
rule fractionalMass:
    input: QCExtractorFractionalMassOutput
    output: fractionalMassOutput
    params: theo_masses=os.path.join(config["R_home"],"theoretical_masses.txt"),r_script=os.path.join(config["R_home"],"fractionalMass.R ")
    run:
        for inp, outp in zip(input,output):
            shell("Rscript {params.r_script} {inp} {outp} {params.theo_masses}")

massAccOutput = createOutputFileNames("massAcc/",QCExtractorIDR38Output,".png")
rule massAcc:
    input: QCExtractorIDR38Output
    output: massAccOutput
    params: r_script=os.path.join(config["R_home"],"MassAcc.R")
    run:
        for inp, outp in zip(input, output):
            shell("Rscript {params.r_script} {inp} {outp}")

massErrOutput = createOutputFileNames("massErr/",QCExtractorIDR38Output,".png")
rule massErr:
    input: QCExtractorIDR38Output
    output: massErrOutput
    params: r_script=os.path.join(config["R_home"],"MassError.R")
    run:
        for inp, outp in zip(input, output):
            shell("Rscript {params.r_script} {inp} {outp}")

ticOutput = createOutputFileNames("tic/",QCExtractorFractionalMassOutput,".png")
rule tic:
    input: QCExtractorTICOutput
    output: ticOutput
    params: r_script=os.path.join(config["R_home"],"TIC.R")
    run:
        for inp, outp in zip(input, output):
            shell("Rscript {params.r_script} {inp} {outp}")


#File with 2 columns is the first one
#second file is with qp QC:0000038 IDR38
#first IDR44, second IDR38
idRatioOutput = createOutputFileNames("idRatio/",QCExtractorIDR38Output,".png")
rule idRatio:
    input: idr38 = QCExtractorIDR38Output, idr44 = QCExtractorIDR44Output
    output: idRatioOutput
    params: r_script=os.path.join(config["R_home"], "IDRatio.R")
    run:
        for idr38, idr44, outp in zip(input.idr38, input.idr44, output):
            shell("Rscript {params.r_script} {idr38} {idr44} {outp}")


QCEmbedderOutput = createOutputFileNames("QCEmbedder/",INPUT_FILES, ".qcML")
rule QCEmbedder:
    input: tic = ticOutput, me = massErrOutput, ma = massAccOutput, fm = fractionalMassOutput,idr = idRatioOutput, _in = QCCalculatorOutput
    output: QCEmbedderOutput
    run:
        for tic, me, ma, fm, idr, _in, outp in zip(input.tic, input.me, input.ma, input.fm, input.idr, input._in, output):
           shell("QCEmbedder -in {_in} -plot {tic} -out {outp} -cv_acc MS:1000235 -qp_att_acc QC:0000023")
           shell("QCEmbedder -in {outp} -plot {me} -out {outp} -cv_acc QC:0000054 -qp_att_acc QC:0000041")
           shell("QCEmbedder -in {outp} -plot {ma} -out {outp} -cv_acc QC:0000053 -qp_att_acc QC:0000041")
           shell("QCEmbedder -in {outp} -plot {fm} -out {outp} -cv_acc QC:0000043 -qp_att_acc QC:0000007")
           shell("QCEmbedder -in {outp} -plot {idr} -out {outp} -cv_acc QC:0000052 -qp_att_acc QC:0000035")


QCShrinkerOutput = RESULT_FILES
rule QCShrinker:
    input: QCEmbedderOutput
    output: QCShrinkerOutput
    run:
        for inp, outp in zip(input, output):
            shell("QCShrinker -in {inp} -out {outp}")
