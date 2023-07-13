a = read.table("/data/yihan/Mxra7_2m/Result_for_X/05.edgeRMarkgenes.Result/T/T_cells.control_treat.glmLRT.plot.txt")
b = read.table("/data/yihan/Mxra7_2m/Result_for_X/05.edgeRMarkgenes.Result/NK/NK_cells.control_treat.glmLRT.plot.txt")
merged_table <- rbind(a, b)
write.table(merged_table, "/data/yihan/Mxra7_2m/Result_for_X/05.edgeRMarkgenes.Result/T_cells.control_treat.glmLRT.plot.txt")
