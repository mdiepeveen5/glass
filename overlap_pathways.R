#needs result_prot_up_gp, result_prot_down_gp, result_up_gp, result_down_gp

pathway_up_both <- result_up_gp %>% 
  inner_join(result_prot_up_gp, by = c("ID"))

#649 RNA, 82 prot, 62 overlap

pathway_up_onlyprot <- result_prot_up_gp %>% 
  anti_join(result_up_gp, by = c("ID"))
