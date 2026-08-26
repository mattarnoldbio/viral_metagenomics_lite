include { ASSIGN_TAXONOMY } from './modules/local/taxonomyassignment/main.nf'

workflow TEST {
    
    main:
    ASSIGN_TAXONOMY (
        params.blast_results,
        params.which_db,
        tuple(params.rankedlineage, params.fullnamelineage)
    )

}