### get MAFs of the selected mutations
odir2 = paste0(odir1,"/MAFs")
system(paste0("mkdir -p ",odir2))

sh_code1 = "/mnt/projects/lailhh/workspace/pipelines/PyClone/PLANET_Pyclone_pipeline/scripts/S01_prepareInputs/S03d_MAF_for_selectedSNPs.sh"

vcfs1 = list.files(odir1,pattern = ".vcf")
vcfs1 = tibble(DNA_lib = gsub(pattern = ".vcf",replacement = "",x = vcfs1),
               vcf = paste0(odir1,vcfs1))
for (v in 1:nrow(vcfs1))
{
  # v = 1
  dna_lib = vcfs1$DNA_lib[v]
  pat = mf %>%
    filter(DNA_lib==dna_lib) %>%
    pull(Unified_ID)

  selectedSNPs = vcfs1$vcf[v]
  pat_tums = mf %>%
    filter(Unified_ID==pat & !is.na(funcotator_dir)) %>%
    pull(DNA_lib)

  lapply(X = pat_tums,
         FUN = function(x) {
           # x = pat_tums[1]
           maf = mf %>%
             filter(DNA_lib==x) %>%
             pull(funcotator_dir)
           cmd = paste0(sh_code1," ",
                        selectedSNPs," ",
                        maf," ",
                        odir2,"/",x,"_selectedSNPs.maf")
           system(cmd)
           return(cmd)
         })

}