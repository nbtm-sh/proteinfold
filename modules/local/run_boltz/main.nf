/*
 * Run Alphafold2
 */
process RUN_BOLTZ {
    tag "$meta.id"
    label 'process_medium'

    container "/srv/scratch/sbf-pipelines/proteinfold/singularity/boltz.sif"
    
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
