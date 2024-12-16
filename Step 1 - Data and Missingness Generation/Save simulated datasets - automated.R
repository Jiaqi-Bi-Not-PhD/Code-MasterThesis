## Save separate data
mar_list <- readRDS("500_mar_lists.RData")
list_names <- names(mar_list)
for(name in list_names) assign(name, mar_list[[name]])
for(name in list_names) saveRDS(mar_list[[name]], file = paste0(name, ".RData"))

