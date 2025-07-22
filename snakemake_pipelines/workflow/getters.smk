def fetch_chrom_str(wildcards, input):
    key = Path(input[0]).stem.split('_')[0]
    if key in list(config['std_chroms']):
        return ' '.join(config['std_chroms'][key])
    else: return ""

def fetch_model(wildcards, input, output):
    ACC = wildcards.acc
    
    models={
        '4kHz': {
            'v4': f"{RERIO_DIR}/res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1",
            'v4_5mCG': f"{DORADO_DIR}/dna_r10.4.1_e8.2_400bps_{ACC}@v4.1.0",
        },
        '5kHz': {
            'v4': f"{DORADO_DIR}/dna_r10.4.1_e8.2_400bps_{ACC}@v4.3.0",
            'v5': f"{DORADO_DIR}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0",
            'v5.2': f"{DORADO_DIR}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.2.0",
        },
    }
    exp, sr, ver, mod = re.search(r'(.*)_([45]kHz)_(sup|hac|fast)_(v[45.2]+)r[123]_?([456]m[CAG]+)?', Path(output[0]).stem).group(1,2,4,5)
    
    if sr=="4" and mod=="5mCG": 
        return models[sr][f"{ver}_{mod}"]
    return models[sr][ver]

def fetch_mod_model(wildcards, input, output):
    ACC = wildcards.acc
    MOD = wildcards.mod
    
    models={
        '4kHz': {
            'v4r2_5mC':  f"{RERIO_DIR}/res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1_5mC@v2",
            'v4r2_6mA':  f"{RERIO_DIR}/res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1_6mA@v2",
            'v4r2_5mCG':  f"{DORADO_DIR}/dna_r10.4.1_e8.2_400bps_{ACC}@v4.1.0_5mCG_5hmCG@v2",
        },
        '5kHz': {
            'v4r1_5mC':  f"{DORADO_DIR}/dna_r10.4.1_e8.2_400bps_{ACC}@v4.3.0_5mC_5hmC@v1",
            'v4r1_5mCG': f"{DORADO_DIR}/dna_r10.4.1_e8.2_400bps_{ACC}@v4.3.0_5mCG_5hmCG@v1",
            'v4r1_6mA':  f"{DORADO_DIR}/dna_r10.4.1_e8.2_400bps_{ACC}@v4.3.0_6mA@v1",
            'v4r1_4mC':  f"{RERIO_DIR}/res_dna_r10.4.1_e8.2_400bps_sup@v4.3.0_4mC_5mC@v1",

            'v4r2_6mA':  f"{DORADO_DIR}/dna_r10.4.1_e8.2_400bps_{ACC}@v4.3.0_6mA@v2",

            'v5r1_5mC':  f"{DORADO_DIR}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0_5mC_5hmC@v1",
            'v5r2_5mC':  f"{DORADO_DIR}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0_5mC_5hmC@v2.0.1",
            'v5r3_5mC':  f"{DORADO_DIR}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0_5mC_5hmC@v3",
            
            'v5r1_5mCG': f"{DORADO_DIR}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0_5mCG_5hmCG@v1",
            'v5r2_5mCG': f"{DORADO_DIR}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0_5mCG_5hmCG@v2.0.1",
            'v5r3_5mCG': f"{DORADO_DIR}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0_5mCG_5hmCG@v3",

            'v5r1_6mA':  f"{DORADO_DIR}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0_6mA@v1",
            'v5r2_6mA':  f"{DORADO_DIR}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0_6mA@v2",
            'v5r3_6mA':  f"{DORADO_DIR}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0_6mA@v3",

            'v5r1_4mC':  f"{DORADO_DIR}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0_4mC_5mC@v1",
            'v5r2_4mC':  f"{DORADO_DIR}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0_4mC_5mC@v2",
            'v5r3_4mC':  f"{DORADO_DIR}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.0.0_4mC_5mC@v3",

            'v5.2r1_5mC':  f"{DORADO_DIR}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.2.0_5mC_5hmC@v1",
            'v5.2r1_5mCG': f"{DORADO_DIR}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.2.0_5mCG_5hmCG@v1",
            'v5.2r1_6mA':  f"{DORADO_DIR}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.2.0_6mA@v1",
            'v5.2r1_4mC':  f"{DORADO_DIR}/dna_r10.4.1_e8.2_400bps_{ACC}@v5.2.0_4mC_5mC@v1",
        }
    }
    exp, sr, ver = re.search(r'(.*)_([45]kHz)_(sup|hac|fast)_(v[45.2]+r[123])', Path(output[0]).stem).group(1,2,4)
    return models[sr][f"{ver}_{MOD}"]

    # return ','.join([ models[sr][f"{ver}_{MOD}"] for MOD in MODS ])

def getRef(wildcards):
    references = config['references']
    return references[wildcards.experiment.split('_')[0]]


def select_dorado(wildcards, input, output):
    sr, acc, ver = re.search(r'.*_([45]kHz)_(sup|hac|fast)_(v[45.2]+r[123])', Path(output[0]).stem).group(1,2,3)
    dorvers = {
        '4kHz_v4r1': '053',
        '4kHz_v4r2': '053',
        '5kHz_v4r1': '091',
        '5kHz_v4r2': '091',
        '5kHz_v5r1': '091',
        '5kHz_v5r2': '091',
        '5kHz_v5r3': '091',
        '5kHz_v5.2r1': '100',
    }
    return dorvers[f'{sr}_{ver}']
