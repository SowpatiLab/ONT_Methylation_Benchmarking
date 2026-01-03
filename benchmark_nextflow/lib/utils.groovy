// general

def fetch_ref(experiment){
    reference_map = new groovy.yaml.YamlSlurper().parseText(file(params.refLookup).text)
    return reference_map['references'][experiment.split('_')[0]]
}

def get_sample_name(bam_file) {
    def fetch_key = "${bam_file}".split('_5kHz')[0]
    return "${fetch_key}_5kHz"
}

def fetch_pod5_loc(bam_file) {
    def fetch_key = get_sample_name(bam_file)
    return "${params.pod5dir}/${fetch_key}"
}

// dorado

def select_base_model(sr, acc, version, mod) {
    def models = [
        '4kHz': [
            'v4': "res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1",
            'v4_5mCG': "dna_r10.4.1_e8.2_400bps_${acc}@v4.1.0",
        ],
        '5kHz': [
            'v4'   : "dna_r10.4.1_e8.2_400bps_${acc}@v4.3.0",
            'v5'   : "dna_r10.4.1_e8.2_400bps_${acc}@v5.0.0",
            'v5.2' : "dna_r10.4.1_e8.2_400bps_${acc}@v5.2.0"
        ]
    ]
    
    if (sr==4 && mod=='5mCG'){
        return models["${sr}kHz_${mod}"]
    } else {
        return models["${sr}kHz"][version.split('r')[0]]
    }
}

def select_mod_model(sr, acc, version, mod){
    def redundant_mods = [
        '_5mC'  : '5mC',
        '_5hmC' : '5mC',
        '_5mCG' : '5mCG',
        '_5hmCG': '5mCG',
        '_4mC'  : '4mC',
        '_6mA'  : '6mA',
    ]

    def models = [
        '4kHz': [
            'v4r2_5mC':  "res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1_5mC@v2",
            'v4r2_6mA':  "res_dna_r10.4.1_e8.2_400bps_sup@v4.0.1_6mA@v2",
            'v4r2_5mCG': "dna_r10.4.1_e8.2_400bps_${acc}@v4.1.0_5mCG_5hmCG@v2",
        ],
        '5kHz': [
            'v4r1_5mC':   "dna_r10.4.1_e8.2_400bps_${acc}@v4.3.0_5mC_5hmC@v1",
            'v5r1_5mC':   "dna_r10.4.1_e8.2_400bps_${acc}@v5.0.0_5mC_5hmC@v1",
            'v5r2_5mC':   "dna_r10.4.1_e8.2_400bps_${acc}@v5.0.0_5mC_5hmC@v2.0.1",
            'v5r3_5mC':   "dna_r10.4.1_e8.2_400bps_${acc}@v5.0.0_5mC_5hmC@v3",
            'v5.2r1_5mC': "dna_r10.4.1_e8.2_400bps_${acc}@v5.2.0_5mC_5hmC@v1",
            'v5.2r2_5mC': "dna_r10.4.1_e8.2_400bps_${acc}@v5.2.0_5mC_5hmC@v2",
            
            'v4r1_5mCG':   "dna_r10.4.1_e8.2_400bps_${acc}@v4.3.0_5mCG_5hmCG@v1",
            'v5r1_5mCG':   "dna_r10.4.1_e8.2_400bps_${acc}@v5.0.0_5mCG_5hmCG@v1",
            'v5r2_5mCG':   "dna_r10.4.1_e8.2_400bps_${acc}@v5.0.0_5mCG_5hmCG@v2.0.1",
            'v5r3_5mCG':   "dna_r10.4.1_e8.2_400bps_${acc}@v5.0.0_5mCG_5hmCG@v3",
            'v5.2r1_5mCG': "dna_r10.4.1_e8.2_400bps_${acc}@v5.2.0_5mCG_5hmCG@v1",
            'v5.2r2_5mCG': "dna_r10.4.1_e8.2_400bps_${acc}@v5.2.0_5mCG_5hmCG@v2",

            'v4r1_6mA':   "dna_r10.4.1_e8.2_400bps_${acc}@v4.3.0_6mA@v1",
            'v4r2_6mA':   "dna_r10.4.1_e8.2_400bps_${acc}@v4.3.0_6mA@v2",
            'v5r1_6mA':   "dna_r10.4.1_e8.2_400bps_${acc}@v5.0.0_6mA@v1",
            'v5r2_6mA':   "dna_r10.4.1_e8.2_400bps_${acc}@v5.0.0_6mA@v2",
            'v5r3_6mA':   "dna_r10.4.1_e8.2_400bps_${acc}@v5.0.0_6mA@v3",
            'v5.2r1_6mA': "dna_r10.4.1_e8.2_400bps_${acc}@v5.2.0_6mA@v1",

            'v4r1_4mC':   "res_dna_r10.4.1_e8.2_400bps_sup@v4.3.0_4mC_5mC@v1",
            'v5r1_4mC':   "dna_r10.4.1_e8.2_400bps_${acc}@v5.0.0_4mC_5mC@v1",
            'v5r2_4mC':   "dna_r10.4.1_e8.2_400bps_${acc}@v5.0.0_4mC_5mC@v2",
            'v5r3_4mC':   "dna_r10.4.1_e8.2_400bps_${acc}@v5.0.0_4mC_5mC@v3",
            'v5.2r1_4mC': "dna_r10.4.1_e8.2_400bps_${acc}@v5.2.0_4mC_5mC@v1",
        ]
    ]
    
    if (mod==null || mod.length()==0) {
        return '--emit-moves'
    } else {
        return models["${sr}kHz"]["${version}_${redundant_mods[mod]}"]
    }
}

// deepmod2

def fetch_deepmod2_model(file, model) {
    models = [
        '5kHz_transformer_v5r3': "transformer_r10.4.1_5khz_v5.0",
        '5kHz_transformer_v4r2': "transformer_r10.4.1_5khz_v4.3",
        '4kHz_transformer_v4r2': "transformer_r10.4.1_4khz_v4.1",
        '5kHz_BiLSTM_v5r3'     : "bilstm_r10.4.1_5khz_v5.0",
        '5kHz_BiLSTM_v4r2'     : "bilstm_r10.4.1_5khz_v4.3",
        '4kHz_BiLSTM_v4r2'     : "bilstm_r10.4.1_4khz_v4.1",
    ]

    pattern = /(.*)_([45]kHz)_(hac|sup)_(v[45\.2]+r[123])/
    matcher = (file =~ pattern)[0]
    sr = matcher[2]
    ver = matcher[4]
    return models["${sr}_${model}_${ver}"]
}

// modkit

def modkitParams(input){
    mod = "${input}".split('.cleansed')[0].split('_')[-1]
    if(mod=="5hmC") { return "--ignore m" }
    else { return " --ignore h"}
}

// samtools

def fetchStdChromosomes(infile) {
    if (!params.filterChromosomes){
        return ""
    } else {
        reference_map = new groovy.yaml.YamlSlurper().parseText(file(params.refLookup).text)
        return reference_map['std_chroms']["${infile}".split('_')[0]].join(' ')
    }
}
