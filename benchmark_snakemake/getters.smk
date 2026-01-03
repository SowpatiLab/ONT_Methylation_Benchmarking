configfile: 'references.yaml'
configfile: 'config.yaml'

Path(f"log").mkdir(parents=True, exist_ok=True)

wildcard_constraints:
    acc='sup|hac|fast',
    mod='(_[564][hm]+[CAG]+|)',
    sr='[45]',
    ver='[45.2]+r[123]',
    deepmodel='transformer|BiLSTM',
    strandstate='(|_stranded)',
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

EXPS = config['exps']
ACC  = config['acc']
RUNS = config['runs']
NRUNS = len(config['runs'])

def fetch_chrom_str(wildcards, input):
    key = Path(input[0]).stem.split('_')[0]
    if key in list(config['std_chroms']):
        return ' '.join(config['std_chroms'][key])
    else: return ""

def fetch_model(wildcards, input, output):
    ACC = wildcards.acc
    
    models={
        '4kHz': {
            'v4': f"{RERIO_MODELS}/res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1",
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
            'v4r2_5mC':  f"{RERIO_MODELS}/res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1_5mC@v2",
            'v4r2_6mA':  f"{RERIO_MODELS}/res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1_6mA@v2",
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

            'v4r1_4mC':   f"{RERIO_MODELS}/res_dna_r10.4.1_e8.2_400bps_sup@v4.3.0_4mC_5mC@v1",
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
    sr, acc, ver = re.search(r'.*_([45]kHz)_(sup|hac|fast)_(v[45.2]+r[123])', Path(output[0]).stem).group(1,2,3)

    dorado_ver = config['dorado_config'][f'{sr}_{ver}'] 
    return config['tooling']['dorado'][dorado_ver]

# tool setup
TOOLS    = config['tooling']
BEDTOOLS = TOOLS['bedtools']
SAMTOOLS = TOOLS['samtools']
MODKIT   = TOOLS['modkit']

# ROCKFISH   = TOOLS['rockfish']
# DEEPMOD2   = TOOLS['deepmod2']
# F5C        = TOOLS['f5c']
# DEEPPLANT  = TOOLS['deepplant']
# DEEPBAM    = TOOLS['deepbam']

#setup dorado version lookup
versions = config['version_model_lookup']
VERSIONS = []
for SR in list(versions):
    for MODEL in versions[SR]:
        VERSIONS += [(SR[0], MODEL[1:])]

RERIO_MODELS  = config['RERIO_MODELS']
DORADO_MODELS = config['DORADO_MODELS']

def fetch_dorado_list():
    if config['dorado_mods']==None: return []
    elif sum([len(config['dorado_mods'][i]) for i in list(config['dorado_mods'])])==0: return []
    
    RUN_MODELS = config['dorado_mods']
    MODS = list(RUN_MODELS)
    RET_LIST = []

    for mod in MODS:
        RET_LIST += expand("meta_context/dorado/{experiment}_5kHz_{acc}_{ver}_{mod}.tsv", experiment=EXPS, acc=ACC, ver=RUN_MODELS[mod], mod=mod)
    RET_LIST = list(filter(lambda x: x.find('hac_v4r1_4mC')<0, RET_LIST))
    return RET_LIST

def fetch_rockfish_list():
    if config['rockfish_run']==None: return []
    return expand("meta_context/rockfish/{experiment}_5kHz_{acc}_{ver}.tsv", experiment=EXPS, acc=ACC, ver=config['rockfish_run'])

def fetch_f5c_list():
    if config['f5c_run']==None: return []
    return expand("meta_context/f5c/{experiment}_5kHz_{acc}_{ver}.tsv", experiment=EXPS, acc=ACC, ver=config['f5c_run'])

def fetch_f5c_stranded_list():
    if config['f5c_stranded_run']==None: return []
    return expand("meta_context/f5c_stranded/{experiment}_5kHz_{acc}_{ver}.tsv", experiment=EXPS, acc=ACC, ver=config['f5c_stranded_run'])

def fetch_deepplant_list():
    if config['deepplnat_run']==None: return []
    return expand("meta_context/deepplant/{experiment}_5kHz_{acc}_{ver}_{context}.tsv", experiment=EXPS, acc=ACC, ver=config['deepplnat_run'], context=['cpg', 'chg'])

def fetch_deepbam_list():
    if config['deepbam_run']==None: return []
    return expand("meta_context/deepbam/{experiment}_5kHz_{acc}_{ver}.tsv", experiment=EXPS, acc=ACC, ver=config['deepbam_run'])


def fetch_deepmod2_bilstm_list():
    if config['deepmod2_bilstm_run']==None: return []
    return expand("meta_context/deepmod2/{experiment}_5kHz_{acc}_{ver}_{deepmodel}.tsv", experiment=EXPS, acc=ACC, deepmodel=['BiLSTM'], ver=config['deepmod2_bilstm_run'])


def fetch_deepmod2_transformer_list():
    if config['deepmod2_transformer_run']==None: return []
    return expand("meta_context/deepmod2/{experiment}_5kHz_{acc}_{ver}_{deepmodel}.tsv", experiment=EXPS, acc=ACC, deepmodel=['transformer'], ver=config['deepmod2_transformer_run'])

# LOAD OTHER TOOLS =================================>
include: 'modules/rockfish.smk'
include: 'modules/deepplant.smk'
include: 'modules/deepbam.smk'
include: 'modules/f5c.smk'
include: 'modules/deepmod2.smk'