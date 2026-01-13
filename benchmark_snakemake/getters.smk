configfile: 'config.json'
configfile: config["refLookup"]

Path(f"log").mkdir(parents=True, exist_ok=True)

wildcard_constraints:
    acc='sup|hac|fast',
    mod='(_[564][hm]+[CAG]+|)',
    sr='[45]',
    ver='[45.2]+r[123]',
    deepmodel='transformer|BiLSTM',
    strandstate='(|_stranded)',
    stranded_ext='(|.stranded)',
    deeptool='deepbam|deepplant',
    modkitver='[45]',
    doradomode='(move|mod)',
    # run=''

timit = bool(config['timit'])
logit = bool(config['logit'])

TIME   = '/usr/bin/time -v -f "%E"' if timit else ''
STDERR = f'2> >(tee -a {{log}} >&2)' if logit else ''

def sh(command):
    final_command = ''.join(command.split('\n'))
    return f"""({TIME} {final_command})  {STDERR}"""

def ntsh(command):
    final_command = ''.join(command.split('\n'))
    return f"""({final_command}) {STDERR}"""

# ==================================================>

EXPS = config['runControls']['experiments']
ACC  = config['runControls']['accuracy']
RUNS = config['runControls']['runs']
NRUNS = len(config['runControls']['runs'])

def fetch_chrom_str(wildcards, input):
    key = Path(input[0]).stem.split('_')[0]
    if key in list(config['std_chroms']):
        return ' '.join(config['std_chroms'][key])
    else: return ""

def fetch_model(wildcards, input, output):
    ACC = wildcards.acc
    
    models={
        '4kHz': {
            'v4': f"{DORADO_MODELS}/res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1",
            'v4_5mCG': f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v4.1.0",
        },
        '5kHz': {
            'v4': f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v4.3.0",
            'v5': f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0",
            'v5.2': f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.2.0",
        },
    }
    exp, sr, ver, mod = re.search(r'(.*)_([45]kHz)_(sup|hac|fast)_(v[45.2]+)r[123]_?([456][hm]+[CAG]+)?', Path(output[0]).stem).group(1,2,4,5)
    
    if sr=="4" and mod=="5mCG": 
        return models[sr][f"{ver}_{mod}"]
    return models[sr][ver]

def fetch_mod_model(wildcards, input, output):
    ACC = wildcards.acc
    if len(wildcards.mod)==0: return ' --emit-moves '
    MOD = {
        '_5mC'  : '5mC',
        '_5hmC' : '5mC',
        '_5mCG' : '5mCG',
        '_5hmCG': '5mCG',
        '_4mC'  : '4mC',
        '_6mA'  : '6mA',
    }[wildcards.mod]

    models={
        '4kHz': {
            'v4r2_5mC':  f"{DORADO_MODELS}/res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1_5mC@v2",
            'v4r2_6mA':  f"{DORADO_MODELS}/res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1_6mA@v2",
            'v4r2_5mCG': f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v4.1.0_5mCG_5hmCG@v2",
        },
        '5kHz': {
            'v4r1_5mC':   f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v4.3.0_5mC_5hmC@v1",
            'v5r1_5mC':   f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0_5mC_5hmC@v1",
            'v5r2_5mC':   f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0_5mC_5hmC@v2.0.1",
            'v5r3_5mC':   f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0_5mC_5hmC@v3",
            'v5.2r1_5mC': f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.2.0_5mC_5hmC@v1",
            'v5.2r2_5mC': f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.2.0_5mC_5hmC@v2",
            
            'v4r1_5mCG':   f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v4.3.0_5mCG_5hmCG@v1",
            'v5r1_5mCG':   f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0_5mCG_5hmCG@v1",
            'v5r2_5mCG':   f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0_5mCG_5hmCG@v2.0.1",
            'v5r3_5mCG':   f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0_5mCG_5hmCG@v3",
            'v5.2r1_5mCG': f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.2.0_5mCG_5hmCG@v1",
            'v5.2r2_5mCG': f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.2.0_5mCG_5hmCG@v2",

            'v4r1_6mA':   f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v4.3.0_6mA@v1",
            'v4r2_6mA':   f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v4.3.0_6mA@v2",
            'v5r1_6mA':   f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0_6mA@v1",
            'v5r2_6mA':   f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0_6mA@v2",
            'v5r3_6mA':   f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0_6mA@v3",
            'v5.2r1_6mA': f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.2.0_6mA@v1",

            'v4r1_4mC':   f"{DORADO_MODELS}/res_dna_r10.4.1_e8.2_400bps_sup@v4.3.0_4mC_5mC@v1",
            'v5r1_4mC':   f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0_4mC_5mC@v1",
            'v5r2_4mC':   f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0_4mC_5mC@v2",
            'v5r3_4mC':   f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0_4mC_5mC@v3",
            'v5.2r1_4mC': f"{DORADO_MODELS}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.2.0_4mC_5mC@v1",
        }
    }
    exp, sr, ver = re.search(r'(.*)_([45]kHz)_(sup|hac|fast)_(v[45.2]+r[123])', Path(output[0]).stem).group(1,2,4)

    return '--modified-bases-models ' + models[sr][f"{ver}_{MOD}"]

def getRef(wildcards):
    references = config['references']
    return references[wildcards.experiment.split('_')[0]]


def select_dorado(wildcards, input, output):
    sr, acc, ver = re.search(r'.*_([45]kHz)_(sup|hac|fast)_(v[45.2]+)r[123]', Path(output[0]).stem).group(1,2,3)

    dorado_ver = config['toolConfig']['dorado']['dorado_version_map'][f'{sr}_{ver}'] 
    return config['toolConfig']['dorado']['dorado_version_aliases'][dorado_ver]

# tool setup
TOOLS    = config['toolConfig']
BEDTOOLS = 'bedtools'
SAMTOOLS = 'samtools'

def fetch_exec(tool, append=True):
    ob = config['toolConfig'][tool]
    
    install_dir = "" if ob.get('install_dir') == None or ob.get('install_dir') == "" else ob.get('install_dir')
    executable  = "" if ob.get('executable') == None or ob.get('executable') == "" else ob.get('executable')
    tooling_dir = "" if config.get('tooling_dir') == None or config.get('tooling_dir') == "" else f"{config.get('tooling_dir')}/"

    if install_dir!="":
        if Path(install_dir).is_absolute():
            return f"{install_dir}/{executable}"
        else:
            if executable == "" or tool=="dorado":
                if Path(f"{tooling_dir}{install_dir}").exists():
                    return f"{tooling_dir}{install_dir}"
            else:
                if Path(f"{tooling_dir}{install_dir}/{executable}").exists():
                    return f"{tooling_dir}{install_dir}/{executable}"
                elif Path(f"{tooling_dir}{install_dir}").exists():
                    return f"{tooling_dir}{install_dir}"
    else:
        return tool

def fetch_model_dir(tool, append=True):
    ob = config['toolConfig'][tool]
    
    model_dir   = config['toolConfig'][tool]['model_dir']
    tooling_dir = "" if config.get('tooling_dir') == None or config.get('tooling_dir') == "" else f"{config.get('tooling_dir')}/"

    if tooling_dir!="":
        return f"{tooling_dir}{model_dir}"
    else:
        return model_dir

toolingDir = "" if config.get('tooling_dir')=="" or config.get('tooling_dir') is None else f"{config['tooling_dir']}/"

DORADO_HOME   = fetch_exec('dorado')
ROCKFISH      = "rockfish"
DEEPMOD2      = fetch_exec('deepmod2')
F5C           = fetch_exec('f5c')
DEEPPLANT     = fetch_exec('deepplant', False)
DEEPBAM       = fetch_exec('deepbam', False)
DORADO_MODELS = fetch_model_dir('dorado')
MODKIT   = fetch_exec('modkit')
ROCKFISH_MODEL = fetch_model_dir('rockfish') + f"/{config['toolConfig']['rockfish']['model']}"

# print(DORADO_HOME)
# print(ROCKFISH)
# print(DEEPMOD2)
# print(F5C)
# print(DEEPPLANT)
# print(DEEPBAM)
# print(DORADO_MODELS)
# print(MODKIT)
# print(DORADO_HOME)

# sys.exit()

#setup dorado version lookup
# versions = config['version_model_lookup']
# VERSIONS = []
# for SR in list(versions):
#     for MODEL in versions[SR]:
#         VERSIONS += [(SR[0], MODEL[1:])]

def fetch_dorado_list():
    if config['runControls']['dorado_mod_models']==None: return []
    elif sum([len(config['runControls']['dorado_mod_models'][i]) for i in list(config['runControls']['dorado_mod_models'])])==0: return []
    
    RUN_MODELS = config['runControls']['dorado_mod_models']
    MODS = list(RUN_MODELS)
    RET_LIST = []

    for mod in MODS:
        RET_LIST += expand("meta/dorado/{experiment}_5kHz_{acc}_{ver}_{mod}.cleansed.ref.std.tsv", experiment=EXPS, acc=ACC, ver=RUN_MODELS[mod], mod=mod)
        # RET_LIST += expand("bam/sorted_mod/{experiment}_5kHz_{acc}_{ver}_{mod}.bam", experiment=EXPS, acc=ACC, ver=RUN_MODELS[mod], mod=mod)
    RET_LIST = list(filter(lambda x: x.find('hac_v4r1_4mC')<0, RET_LIST))
    return RET_LIST

def fetch_rockfish_list():
    if config['runControls']['other_tools']['rockfish']==None: return []
    return expand("meta/rockfish/{experiment}_5kHz_{acc}_{ver}.rockfish.aggregated.rebed.ref.std.bed", experiment=EXPS, acc=ACC, ver=config['runControls']['other_tools']['rockfish'])

def fetch_f5c_list():
    if config['runControls']['other_tools']['f5c']==None: return []
    return expand("meta/f5c/{experiment}_5kHz_{acc}_{ver}.f5c.aggregated.rebed.ref.std.bed", experiment=EXPS, acc=ACC, ver=config['runControls']['other_tools']['f5c'])

def fetch_f5c_stranded_list():
    if config['runControls']['other_tools']['f5c_stranded']==None: return []
    return expand("meta/f5c_stranded/{experiment}_5kHz_{acc}_{ver}.f5c.stranded.aggregated.rebed.ref.std.bed", experiment=EXPS, acc=ACC, ver=config['runControls']['other_tools']['f5c_stranded'])

def fetch_deepplant_list():
    if config['runControls']['other_tools']['deepplant']==None: return []
    return expand("meta/deepplant/{experiment}_5kHz_{acc}_{ver}.deepplant_{context}.aggregated.rebed.ref.std.bed", experiment=EXPS, acc=ACC, ver=config['runControls']['other_tools']['deepplant'], context=['cpg', 'chg'])

def fetch_deepbam_list():
    if config['runControls']['other_tools']['deepbam']==None: return []
    return expand("meta/deepbam/{experiment}_5kHz_{acc}_{ver}.deepbam.aggregated.rebed.ref.std.bed", experiment=EXPS, acc=ACC, ver=config['runControls']['other_tools']['deepbam'])


def fetch_deepmod2_bilstm_list():
    if config['runControls']['other_tools']['deepmod2_bilstm']==None: return []
    return expand("meta/deepmod2/{experiment}_5kHz_{acc}_{ver}_{deepmodel}.deepmod2.aggregated.rebed.ref.std.bed", experiment=EXPS, acc=ACC, deepmodel=['BiLSTM'], ver=config['runControls']['other_tools']['deepmod2_bilstm'])


def fetch_deepmod2_transformer_list():
    if config['runControls']['other_tools']['deepmod2_transformer']==None: return []
    return expand("meta/deepmod2/{experiment}_5kHz_{acc}_{ver}_{deepmodel}.deepmod2.aggregated.rebed.ref.std.bed", experiment=EXPS, acc=ACC, deepmodel=['transformer'], ver=config['runControls']['other_tools']['deepmod2_transformer'])


## QC Controls 
def select_QC_default_model():
    models = []
    other_tools = config['runControls']['other_tools']
    for tool in list(other_tools):
        models += other_tools[tool]
    models = list(set(models))

    if len(models)==0: return [ config['runControls']['qc']['default_model'] ]
    else: return models

def fetch_qc_tool_list(tool):
    ext = {
        "nanostat": ".nanostat",
        "nanoplot": "",
        "mosdepth": "",
        "nanoq": ".nanoq"
    }
    if config['runControls']['qc']['tools'][tool]: return []
    else: return expand(f"qc/{tool}/{{experiment}}_5kHz_{{acc}}_{{ver}}{ext[tool]}", experiment=EXPS, acc=ACC, ver=select_QC_default_model())


# LOAD OTHER TOOLS =================================>
include: 'modules_smk/rockfish.smk'
include: 'modules_smk/deepplant.smk'
include: 'modules_smk/deepbam.smk'
include: 'modules_smk/f5c.smk'
include: 'modules_smk/deepmod2.smk'
