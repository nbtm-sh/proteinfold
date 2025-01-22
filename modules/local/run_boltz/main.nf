/*
 * Run Alphafold2
 */
process RUN_BOLTZ {
    tag "$meta.id"
    label 'process_medium'

    container "/home/z3545907/boltz_apptainer/boltz.sif"
    
    input:
    tuple val(meta), path(fasta)
    path ('boltz1_conf.ckpt')
    path ('ccd.pkl')
    
    output:
    path ("boltz_results_${fasta.baseName}/processed/msa/*.npz"), emit: msa
    path ("boltz_results_${fasta.baseName}/processed/structures/*.npz"), emit: structures
    path ("boltz_results_${fasta.baseName}/predictions/${fasta.baseName}/confidence*.json"), emit: confidence
    path ("boltz_results_${fasta.baseName}/predictions/${fasta.baseName}/plddt_*.npz"), emit: plddt
    
    script:
    """
    /opt/miniforge/envs/boltz/bin/boltz predict --use_msa_server --accelerator gpu "./${fasta.name}" --cache ./
    """
}

process RUN_ALPHAFOLD2 {
    tag "$meta.id"
    label 'process_medium'

    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error("Local RUN_ALPHAFOLD2 module does not support Conda. Please use Docker / Singularity / Podman instead.")
    }

    container "nf-core/proteinfold_alphafold2_standard:dev"

    input:
    tuple val(meta), path(fasta)
    val   db_preset
    val   alphafold2_model_preset
    path ('params/*')
    path ('bfd/*')
    path ('small_bfd/*')
    path ('mgnify/*')
    path ('pdb70/*')
    path ('pdb_mmcif/*')
    path ('uniref30/*')
    path ('uniref90/*')
    path ('pdb_seqres/*')
    path ('uniprot/*')

    output:
    path ("${fasta.baseName}*")
    tuple val(meta), path ("${fasta.baseName}/ranked*pdb"), emit: pdb
    tuple val(meta), path ("${fasta.baseName}/*_msa.tsv") , emit: msa
    tuple val(meta), path ("*_mqc.tsv")                   , emit: multiqc
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def db_preset = db_preset ? "full_dbs --bfd_database_path=${params.alphafold2_db}/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt --uniref30_database_path=${params.alphafold2_db}/uniref30/UniRef30_2021_03" :
        "reduced_dbs --small_bfd_database_path=${params.alphafold2_db}/small_bfd/bfd-first_non_consensus_sequences.fasta"
    if (alphafold2_model_preset == 'multimer') {
        alphafold2_model_preset += " --pdb_seqres_database_path=${params.alphafold2_db}/pdb_seqres/pdb_seqres.txt --uniprot_database_path=${params.alphafold2_db}/uniprot/uniprot.fasta "
    }
    else {
        alphafold2_model_preset += " --pdb70_database_path=${params.alphafold2_db}/pdb70/pdb70_from_mmcif_200916/pdb70 "
    }
    """
    if [ -d ${params.alphafold2_db}/params/ ]; then ln -r -s ${params.alphafold2_db}/params params; fi
    python3 /app/alphafold/run_alphafold.py \
        --fasta_paths=${fasta} \
        --model_preset=${alphafold2_model_preset} \
        --db_preset=${db_preset} \
        --output_dir=\$PWD \
        --data_dir=\$PWD \
        --uniref90_database_path=${params.alphafold2_db}/uniref90/uniref90.fasta \
        --mgnify_database_path=${params.alphafold2_db}/mgnify/mgy_clusters_2022_05.fa \
        --template_mmcif_dir=${params.alphafold2_db}/pdb_mmcif/mmcif_files \
        --obsolete_pdbs_path=${params.alphafold2_db}/pdb_mmcif/obsolete.dat \
        --use_gpu_relax \
        $args

    cp "${fasta.baseName}"/ranked_0.pdb ./"${fasta.baseName}".alphafold.pdb
    cd "${fasta.baseName}"
    awk '{print \$6"\\t"\$11}' ranked_0.pdb | uniq > ranked_0_plddt.tsv
    for i in 1 2 3 4
        do awk '{print \$6"\\t"\$11}' ranked_\$i.pdb | uniq | awk '{print \$2}' > ranked_"\$i"_plddt.tsv
    done
    paste ranked_0_plddt.tsv ranked_1_plddt.tsv ranked_2_plddt.tsv ranked_3_plddt.tsv ranked_4_plddt.tsv > plddt.tsv
    echo -e Positions"\\t"rank_0"\\t"rank_1"\\t"rank_2"\\t"rank_3"\\t"rank_4 > header.tsv
    cat header.tsv plddt.tsv > ../"${fasta.baseName}"_plddt_mqc.tsv

    extract_output.py --name ${fasta.baseName} \\
        --pkls features.pkl
    cd ..

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch ./"${fasta.baseName}".alphafold.pdb
    touch ./"${fasta.baseName}"_mqc.tsv
    mkdir "${fasta.baseName}"
    touch "${fasta.baseName}/ranked_0.pdb"
    touch "${fasta.baseName}/ranked_1.pdb"
    touch "${fasta.baseName}/ranked_2.pdb"
    touch "${fasta.baseName}/ranked_3.pdb"
    touch "${fasta.baseName}/ranked_4.pdb"
    touch "${fasta.baseName}/${fasta.baseName}_msa.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
