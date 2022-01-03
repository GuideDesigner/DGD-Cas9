library(data.table)
data_connection = read.csv("CC_feature.csv",header = TRUE)
#head(data_connection)

size_1 = subset(data_connection,data_connection[,ncol(data_connection)]==1)
if (length(row.names(size_1)) > 0)
{
  serial_id_list = unique(size_1[,4])
  serial_id_1 = NULL
  loop = NULL
  for (i in 1:length(serial_id_list))
  {
    serial_id_1[[i]] = subset(size_1,size_1[,4]==serial_id_list[[i]])
    serial_id_1[[i]]["Unique_ID"] = paste0(serial_id_1[[i]][1,2],"_",serial_id_1[[i]][1,3],"_",serial_id_1[[i]][1,4],sep="")
  }
  size_1_data = as.data.frame(rbindlist(serial_id_1,fill=TRUE))
  
  for (i in 1:length(rownames(size_1_data)))
  {
    if (size_1_data[i,3]>=32 & size_1_data[i,3]<=37)
    {
      loop[[i]] = "TL"
    }
    else if (size_1_data[i,3]>=54 & size_1_data[i,3]<=57)
    {
      loop[[i]] = "SL1"
    }
    else if (size_1_data[i,3]>=72 & size_1_data[i,3]<=77)
    {
      loop[[i]] = "SL2"
    }
    else if (size_1_data[i,3]>=88 & size_1_data[i,3]<=90)
    {
      loop[[i]] = "SL3"
    }
    else if (size_1_data[i,3]>=21 & size_1_data[i,3]<=31)
    {
      loop[[i]] = "R"
    }
    else if (size_1_data[i,3]>=38 & size_1_data[i,3]<=53)
    {
      loop[[i]] = "AR"
    }
    else if (size_1_data[i,3]>=63 & size_1_data[i,3]<=67)
    {
      loop[[i]] = "LR"
    }
    else
    {
      loop[[i]] = "NS"
    }
  }
  size_1_data["Structure"] = loop  
}

size_2 = subset(data_connection,data_connection[,ncol(data_connection)]==2)
if (length(row.names(size_2)) > 0)
{
  serial_id_list = unique(size_2[,4])
  serial_id_2 = NULL
  loop = NULL
  for (i in 1:length(serial_id_list))
  {
    serial_id_2[[i]] = subset(size_2,size_2[,4]==serial_id_list[[i]])
    serial_id_2[[i]]["Unique_ID"] = paste0(serial_id_2[[i]][1,2],"_",serial_id_2[[i]][1,3],"_",serial_id_2[[i]][1,4],sep="")
  }
  size_2_data = as.data.frame(rbindlist(serial_id_2,fill=TRUE))
  
  for (i in 1:length(rownames(size_2_data)))
  {
    if (size_2_data[i,3]>=32 & size_2_data[i,3]<=37)
    {
      loop[[i]] = "TL"
    }
    else if (size_2_data[i,3]>=54 & size_2_data[i,3]<=57)
    {
      loop[[i]] = "SL1"
    }
    else if (size_2_data[i,3]>=72 & size_2_data[i,3]<=77)
    {
      loop[[i]] = "SL2"
    }
    else if (size_2_data[i,3]>=88 & size_2_data[i,3]<=90)
    {
      loop[[i]] = "SL3"
    }
    else if (size_2_data[i,3]>=21 & size_2_data[i,3]<=31)
    {
      loop[[i]] = "R"
    }
    else if (size_2_data[i,3]>=38 & size_2_data[i,3]<=53)
    {
      loop[[i]] = "AR"
    }
    else if (size_2_data[i,3]>=63 & size_2_data[i,3]<=67)
    {
      loop[[i]] = "LR"
    }
    else
    {
      loop[[i]] = "NS"
    }
  }
  size_2_data["Structure"] = loop  
}

size_3 = subset(data_connection,data_connection[,ncol(data_connection)]==3)
if (length(row.names(size_3)) > 0)
{
  serial_id_list = unique(size_3[,4])
  serial_id_3 = NULL
  loop = NULL
  for (i in 1:length(serial_id_list))
  {
    serial_id_3[[i]] = subset(size_3,size_3[,4]==serial_id_list[[i]])
    serial_id_3[[i]]["Unique_ID"] = paste0(serial_id_3[[i]][1,2],"_",serial_id_3[[i]][1,3],"_",serial_id_3[[i]][1,4],sep="")
  }
  size_3_data = as.data.frame(rbindlist(serial_id_3,fill=TRUE))
  
  for (i in 1:length(rownames(size_3_data)))
  {
    if (size_3_data[i,3]>=32 & size_3_data[i,3]<=37)
    {
      loop[[i]] = "TL"
    }
    else if (size_3_data[i,3]>=54 & size_3_data[i,3]<=57)
    {
      loop[[i]] = "SL1"
    }
    else if (size_3_data[i,3]>=72 & size_3_data[i,3]<=77)
    {
      loop[[i]] = "SL2"
    }
    else if (size_3_data[i,3]>=88 & size_3_data[i,3]<=90)
    {
      loop[[i]] = "SL3"
    }
    else if (size_3_data[i,3]>=21 & size_3_data[i,3]<=31)
    {
      loop[[i]] = "R"
    }
    else if (size_3_data[i,3]>=38 & size_3_data[i,3]<=53)
    {
      loop[[i]] = "AR"
    }
    else if (size_3_data[i,3]>=63 & size_3_data[i,3]<=67)
    {
      loop[[i]] = "LR"
    }
    else
    {
      loop[[i]] = "NS"
    }
  }
  size_3_data["Structure"] = loop  
}

size_4 = subset(data_connection,data_connection[,ncol(data_connection)]==4)
if (length(row.names(size_4)) > 0)
{
  serial_id_list = unique(size_4[,4])
  serial_id_4 = NULL
  loop = NULL
  for (i in 1:length(serial_id_list))
  {
    serial_id_4[[i]] = subset(size_4,size_4[,4]==serial_id_list[[i]])
    serial_id_4[[i]]["Unique_ID"] = paste0(serial_id_4[[i]][1,2],"_",serial_id_4[[i]][1,3],"_",serial_id_4[[i]][1,4],sep="")
  }
  size_4_data = as.data.frame(rbindlist(serial_id_4,fill=TRUE))
  
  for (i in 1:length(rownames(size_4_data)))
  {
    if (size_4_data[i,3]>=32 & size_4_data[i,3]<=37)
    {
      loop[[i]] = "TL"
    }
    else if (size_4_data[i,3]>=54 & size_4_data[i,3]<=57)
    {
      loop[[i]] = "SL1"
    }
    else if (size_4_data[i,3]>=72 & size_4_data[i,3]<=77)
    {
      loop[[i]] = "SL2"
    }
    else if (size_4_data[i,3]>=88 & size_4_data[i,3]<=90)
    {
      loop[[i]] = "SL3"
    }
    else if (size_4_data[i,3]>=21 & size_4_data[i,3]<=31)
    {
      loop[[i]] = "R"
    }
    else if (size_4_data[i,3]>=38 & size_4_data[i,3]<=53)
    {
      loop[[i]] = "AR"
    }
    else if (size_4_data[i,3]>=63 & size_4_data[i,3]<=67)
    {
      loop[[i]] = "LR"
    }
    else
    {
      loop[[i]] = "NS"
    }
  }
  size_4_data["Structure"] = loop  
}

size_5 = subset(data_connection,data_connection[,ncol(data_connection)]==5)
if (length(row.names(size_5)) > 0)
{
  serial_id_list = unique(size_5[,4])
  serial_id_5 = NULL
  loop = NULL
  for (i in 1:length(serial_id_list))
  {
    serial_id_5[[i]] = subset(size_5,size_5[,4]==serial_id_list[[i]])
    serial_id_5[[i]]["Unique_ID"] = paste0(serial_id_5[[i]][1,2],"_",serial_id_5[[i]][1,3],"_",serial_id_5[[i]][1,4],sep="")
  }
  size_5_data = as.data.frame(rbindlist(serial_id_5,fill=TRUE))
  
  for (i in 1:length(rownames(size_5_data)))
  {
    if (size_5_data[i,3]>=32 & size_5_data[i,3]<=37)
    {
      loop[[i]] = "TL"
    }
    else if (size_5_data[i,3]>=54 & size_5_data[i,3]<=57)
    {
      loop[[i]] = "SL1"
    }
    else if (size_5_data[i,3]>=72 & size_5_data[i,3]<=77)
    {
      loop[[i]] = "SL2"
    }
    else if (size_5_data[i,3]>=88 & size_5_data[i,3]<=90)
    {
      loop[[i]] = "SL3"
    }
    else if (size_5_data[i,3]>=21 & size_5_data[i,3]<=31)
    {
      loop[[i]] = "R"
    }
    else if (size_5_data[i,3]>=38 & size_5_data[i,3]<=53)
    {
      loop[[i]] = "AR"
    }
    else if (size_5_data[i,3]>=63 & size_5_data[i,3]<=67)
    {
      loop[[i]] = "LR"
    }
    else
    {
      loop[[i]] = "NS"
    }
  }
  size_5_data["Structure"] = loop  
}

size_6 = subset(data_connection,data_connection[,ncol(data_connection)]==6)
if (length(row.names(size_6)) > 0)
{
  serial_id_list = unique(size_6[,4])
  serial_id_6 = NULL
  loop = NULL
  for (i in 1:length(serial_id_list))
  {
    serial_id_6[[i]] = subset(size_6,size_6[,4]==serial_id_list[[i]])
    serial_id_6[[i]]["Unique_ID"] = paste0(serial_id_6[[i]][1,2],"_",serial_id_6[[i]][1,3],"_",serial_id_6[[i]][1,4],sep="")
  }
  size_6_data = as.data.frame(rbindlist(serial_id_6,fill=TRUE))
  
  for (i in 1:length(rownames(size_6_data)))
  {
    if (size_6_data[i,3]>=32 & size_6_data[i,3]<=37)
    {
      loop[[i]] = "TL"
    }
    else if (size_6_data[i,3]>=54 & size_6_data[i,3]<=57)
    {
      loop[[i]] = "SL1"
    }
    else if (size_6_data[i,3]>=72 & size_6_data[i,3]<=77)
    {
      loop[[i]] = "SL2"
    }
    else if (size_6_data[i,3]>=88 & size_6_data[i,3]<=90)
    {
      loop[[i]] = "SL3"
    }
    else if (size_6_data[i,3]>=21 & size_6_data[i,3]<=31)
    {
      loop[[i]] = "R"
    }
    else if (size_6_data[i,3]>=38 & size_6_data[i,3]<=53)
    {
      loop[[i]] = "AR"
    }
    else if (size_6_data[i,3]>=63 & size_6_data[i,3]<=67)
    {
      loop[[i]] = "LR"
    }
    else
    {
      loop[[i]] = "NS"
    }
  }
  size_6_data["Structure"] = loop  
}

size_7 = subset(data_connection,data_connection[,ncol(data_connection)]==7)
if (length(row.names(size_7)) > 0)
{
  serial_id_list = unique(size_7[,4])
  serial_id_7 = NULL
  loop = NULL
  for (i in 1:length(serial_id_list))
  {
    serial_id_7[[i]] = subset(size_7,size_7[,4]==serial_id_list[[i]])
    serial_id_7[[i]]["Unique_ID"] = paste0(serial_id_7[[i]][1,2],"_",serial_id_7[[i]][1,3],"_",serial_id_7[[i]][1,4],sep="")
  }
  size_7_data = as.data.frame(rbindlist(serial_id_7,fill=TRUE))
  
  for (i in 1:length(rownames(size_7_data)))
  {
    if (size_7_data[i,3]>=32 & size_7_data[i,3]<=37)
    {
      loop[[i]] = "TL"
    }
    else if (size_7_data[i,3]>=54 & size_7_data[i,3]<=57)
    {
      loop[[i]] = "SL1"
    }
    else if (size_7_data[i,3]>=72 & size_7_data[i,3]<=77)
    {
      loop[[i]] = "SL2"
    }
    else if (size_7_data[i,3]>=88 & size_7_data[i,3]<=90)
    {
      loop[[i]] = "SL3"
    }
    else if (size_7_data[i,3]>=21 & size_7_data[i,3]<=31)
    {
      loop[[i]] = "R"
    }
    else if (size_7_data[i,3]>=38 & size_7_data[i,3]<=53)
    {
      loop[[i]] = "AR"
    }
    else if (size_7_data[i,3]>=63 & size_7_data[i,3]<=67)
    {
      loop[[i]] = "LR"
    }
    else
    {
      loop[[i]] = "NS"
    }
  }
  size_7_data["Structure"] = loop  
}

size_8 = subset(data_connection,data_connection[,ncol(data_connection)]==8)
if (length(row.names(size_8)) > 0)
{
  serial_id_list = unique(size_8[,4])
  serial_id_8 = NULL
  loop = NULL
  for (i in 1:length(serial_id_list))
  {
    serial_id_8[[i]] = subset(size_8,size_8[,4]==serial_id_list[[i]])
    serial_id_8[[i]]["Unique_ID"] = paste0(serial_id_8[[i]][1,2],"_",serial_id_8[[i]][1,3],"_",serial_id_8[[i]][1,4],sep="")
  }
  size_8_data = as.data.frame(rbindlist(serial_id_8,fill=TRUE))
  
  for (i in 1:length(rownames(size_8_data)))
  {
    if (size_8_data[i,3]>=32 & size_8_data[i,3]<=37)
    {
      loop[[i]] = "TL"
    }
    else if (size_8_data[i,3]>=54 & size_8_data[i,3]<=57)
    {
      loop[[i]] = "SL1"
    }
    else if (size_8_data[i,3]>=72 & size_8_data[i,3]<=77)
    {
      loop[[i]] = "SL2"
    }
    else if (size_8_data[i,3]>=88 & size_8_data[i,3]<=90)
    {
      loop[[i]] = "SL3"
    }
    else if (size_8_data[i,3]>=21 & size_8_data[i,3]<=31)
    {
      loop[[i]] = "R"
    }
    else if (size_8_data[i,3]>=38 & size_8_data[i,3]<=53)
    {
      loop[[i]] = "AR"
    }
    else if (size_8_data[i,3]>=63 & size_8_data[i,3]<=67)
    {
      loop[[i]] = "LR"
    }
    else
    {
      loop[[i]] = "NS"
    }
  }
  size_8_data["Structure"] = loop  
}

size_9 = subset(data_connection,data_connection[,ncol(data_connection)]==9)
if (length(row.names(size_9)) > 0)
{
  serial_id_list = unique(size_9[,4])
  serial_id_9 = NULL
  loop = NULL
  for (i in 1:length(serial_id_list))
  {
    serial_id_9[[i]] = subset(size_9,size_9[,4]==serial_id_list[[i]])
    serial_id_9[[i]]["Unique_ID"] = paste0(serial_id_9[[i]][1,2],"_",serial_id_9[[i]][1,3],"_",serial_id_9[[i]][1,4],sep="")
  }
  size_9_data = as.data.frame(rbindlist(serial_id_9,fill=TRUE))
  
  for (i in 1:length(rownames(size_9_data)))
  {
    if (size_9_data[i,3]>=32 & size_9_data[i,3]<=37)
    {
      loop[[i]] = "TL"
    }
    else if (size_9_data[i,3]>=54 & size_9_data[i,3]<=57)
    {
      loop[[i]] = "SL1"
    }
    else if (size_9_data[i,3]>=72 & size_9_data[i,3]<=77)
    {
      loop[[i]] = "SL2"
    }
    else if (size_9_data[i,3]>=88 & size_9_data[i,3]<=90)
    {
      loop[[i]] = "SL3"
    }
    else if (size_9_data[i,3]>=21 & size_9_data[i,3]<=31)
    {
      loop[[i]] = "R"
    }
    else if (size_9_data[i,3]>=38 & size_9_data[i,3]<=53)
    {
      loop[[i]] = "AR"
    }
    else if (size_9_data[i,3]>=63 & size_9_data[i,3]<=67)
    {
      loop[[i]] = "LR"
    }
    else
    {
      loop[[i]] = "NS"
    }
  }
  size_9_data["Structure"] = loop  
}

size_10 = subset(data_connection,data_connection[,ncol(data_connection)]==10)
if (length(row.names(size_10)) > 0)
{
  serial_id_list = unique(size_10[,4])
  serial_id_10 = NULL
  loop = NULL
  for (i in 1:length(serial_id_list))
  {
    serial_id_10[[i]] = subset(size_10,size_10[,4]==serial_id_list[[i]])
    serial_id_10[[i]]["Unique_ID"] = paste0(serial_id_10[[i]][1,2],"_",serial_id_10[[i]][1,3],"_",serial_id_10[[i]][1,4],sep="")
  }
  size_10_data = as.data.frame(rbindlist(serial_id_10,fill=TRUE))
  
  for (i in 1:length(rownames(size_10_data)))
  {
    if (size_10_data[i,3]>=32 & size_10_data[i,3]<=37)
    {
      loop[[i]] = "TL"
    }
    else if (size_10_data[i,3]>=54 & size_10_data[i,3]<=57)
    {
      loop[[i]] = "SL1"
    }
    else if (size_10_data[i,3]>=72 & size_10_data[i,3]<=77)
    {
      loop[[i]] = "SL2"
    }
    else if (size_10_data[i,3]>=88 & size_10_data[i,3]<=90)
    {
      loop[[i]] = "SL3"
    }
    else if (size_10_data[i,3]>=21 & size_10_data[i,3]<=31)
    {
      loop[[i]] = "R"
    }
    else if (size_10_data[i,3]>=38 & size_10_data[i,3]<=53)
    {
      loop[[i]] = "AR"
    }
    else if (size_10_data[i,3]>=63 & size_10_data[i,3]<=67)
    {
      loop[[i]] = "LR"
    }
    else
    {
      loop[[i]] = "NS"
    }
  }
  size_10_data["Structure"] = loop  
}

size_11 = subset(data_connection,data_connection[,ncol(data_connection)]==11)
if (length(row.names(size_11)) > 0)
{
  serial_id_list = unique(size_11[,4])
  serial_id_11 = NULL
  loop = NULL
  for (i in 1:length(serial_id_list))
  {
    serial_id_11[[i]] = subset(size_11,size_11[,4]==serial_id_list[[i]])
    serial_id_11[[i]]["Unique_ID"] = paste0(serial_id_11[[i]][1,2],"_",serial_id_11[[i]][1,3],"_",serial_id_11[[i]][1,4],sep="")
  }
  size_11_data = as.data.frame(rbindlist(serial_id_11,fill=TRUE))
  
  for (i in 1:length(rownames(size_11_data)))
  {
    if (size_11_data[i,3]>=32 & size_11_data[i,3]<=37)
    {
      loop[[i]] = "TL"
    }
    else if (size_11_data[i,3]>=54 & size_11_data[i,3]<=57)
    {
      loop[[i]] = "SL1"
    }
    else if (size_11_data[i,3]>=72 & size_11_data[i,3]<=77)
    {
      loop[[i]] = "SL2"
    }
    else if (size_11_data[i,3]>=88 & size_11_data[i,3]<=90)
    {
      loop[[i]] = "SL3"
    }
    else if (size_11_data[i,3]>=21 & size_11_data[i,3]<=31)
    {
      loop[[i]] = "R"
    }
    else if (size_11_data[i,3]>=38 & size_11_data[i,3]<=53)
    {
      loop[[i]] = "AR"
    }
    else if (size_11_data[i,3]>=63 & size_11_data[i,3]<=67)
    {
      loop[[i]] = "LR"
    }
    else
    {
      loop[[i]] = "NS"
    }
  }
  size_11_data["Structure"] = loop  
}
size_12 = subset(data_connection,data_connection[,ncol(data_connection)]==12)
if (length(row.names(size_12)) > 0)
{
  serial_id_list = unique(size_12[,4])
  serial_id_12 = NULL
  loop = NULL
  for (i in 1:length(serial_id_list))
  {
    serial_id_12[[i]] = subset(size_12,size_12[,4]==serial_id_list[[i]])
    serial_id_12[[i]]["Unique_ID"] = paste0(serial_id_12[[i]][1,2],"_",serial_id_12[[i]][1,3],"_",serial_id_12[[i]][1,4],sep="")
  }
  size_12_data = as.data.frame(rbindlist(serial_id_12,fill=TRUE))
  
  for (i in 1:length(rownames(size_12_data)))
  {
    if (size_12_data[i,3]>=32 & size_12_data[i,3]<=37)
    {
      loop[[i]] = "TL"
    }
    else if (size_12_data[i,3]>=54 & size_12_data[i,3]<=57)
    {
      loop[[i]] = "SL1"
    }
    else if (size_12_data[i,3]>=72 & size_12_data[i,3]<=77)
    {
      loop[[i]] = "SL2"
    }
    else if (size_12_data[i,3]>=88 & size_12_data[i,3]<=90)
    {
      loop[[i]] = "SL3"
    }
    else if (size_12_data[i,3]>=21 & size_12_data[i,3]<=31)
    {
      loop[[i]] = "R"
    }
    else if (size_12_data[i,3]>=38 & size_12_data[i,3]<=53)
    {
      loop[[i]] = "AR"
    }
    else if (size_12_data[i,3]>=63 & size_12_data[i,3]<=67)
    {
      loop[[i]] = "LR"
    }
    else
    {
      loop[[i]] = "NS"
    }
  }
  size_12_data["Structure"] = loop  
}
size_13 = subset(data_connection,data_connection[,ncol(data_connection)]==13)
if (length(row.names(size_13)) > 0)
{
  serial_id_list = unique(size_13[,4])
  serial_id_13 = NULL
  loop = NULL
  for (i in 1:length(serial_id_list))
  {
    serial_id_13[[i]] = subset(size_13,size_13[,4]==serial_id_list[[i]])
    serial_id_13[[i]]["Unique_ID"] = paste0(serial_id_13[[i]][1,2],"_",serial_id_13[[i]][1,3],"_",serial_id_13[[i]][1,4],sep="")
  }
  size_13_data = as.data.frame(rbindlist(serial_id_13,fill=TRUE))
  
  for (i in 1:length(rownames(size_13_data)))
  {
    if (size_13_data[i,3]>=32 & size_13_data[i,3]<=37)
    {
      loop[[i]] = "TL"
    }
    else if (size_13_data[i,3]>=54 & size_13_data[i,3]<=57)
    {
      loop[[i]] = "SL1"
    }
    else if (size_13_data[i,3]>=72 & size_13_data[i,3]<=77)
    {
      loop[[i]] = "SL2"
    }
    else if (size_13_data[i,3]>=88 & size_13_data[i,3]<=90)
    {
      loop[[i]] = "SL3"
    }
    else if (size_13_data[i,3]>=21 & size_13_data[i,3]<=31)
    {
      loop[[i]] = "R"
    }
    else if (size_13_data[i,3]>=38 & size_13_data[i,3]<=53)
    {
      loop[[i]] = "AR"
    }
    else if (size_13_data[i,3]>=63 & size_13_data[i,3]<=67)
    {
      loop[[i]] = "LR"
    }
    else
    {
      loop[[i]] = "NS"
    }
  }
  size_13_data["Structure"] = loop  
}
size_14 = subset(data_connection,data_connection[,ncol(data_connection)]==14)
if (length(row.names(size_14)) > 0)
{
  serial_id_list = unique(size_14[,4])
  serial_id_14 = NULL
  loop = NULL
  for (i in 1:length(serial_id_list))
  {
    serial_id_14[[i]] = subset(size_14,size_14[,4]==serial_id_list[[i]])
    serial_id_14[[i]]["Unique_ID"] = paste0(serial_id_14[[i]][1,2],"_",serial_id_14[[i]][1,3],"_",serial_id_14[[i]][1,4],sep="")
  }
  size_14_data = as.data.frame(rbindlist(serial_id_14,fill=TRUE))
  
  for (i in 1:length(rownames(size_14_data)))
  {
    if (size_14_data[i,3]>=32 & size_14_data[i,3]<=37)
    {
      loop[[i]] = "TL"
    }
    else if (size_14_data[i,3]>=54 & size_14_data[i,3]<=57)
    {
      loop[[i]] = "SL1"
    }
    else if (size_14_data[i,3]>=72 & size_14_data[i,3]<=77)
    {
      loop[[i]] = "SL2"
    }
    else if (size_14_data[i,3]>=88 & size_14_data[i,3]<=90)
    {
      loop[[i]] = "SL3"
    }
    else if (size_14_data[i,3]>=21 & size_14_data[i,3]<=31)
    {
      loop[[i]] = "R"
    }
    else if (size_14_data[i,3]>=38 & size_14_data[i,3]<=53)
    {
      loop[[i]] = "AR"
    }
    else if (size_14_data[i,3]>=63 & size_14_data[i,3]<=67)
    {
      loop[[i]] = "LR"
    }
    else
    {
      loop[[i]] = "NS"
    }
  }
  size_14_data["Structure"] = loop  
}
size_15 = subset(data_connection,data_connection[,ncol(data_connection)]==15)
if (length(row.names(size_15)) > 0)
{
  serial_id_list = unique(size_15[,4])
  serial_id_15 = NULL
  loop = NULL
  for (i in 1:length(serial_id_list))
  {
    serial_id_15[[i]] = subset(size_15,size_15[,4]==serial_id_list[[i]])
    serial_id_15[[i]]["Unique_ID"] = paste0(serial_id_15[[i]][1,2],"_",serial_id_15[[i]][1,3],"_",serial_id_15[[i]][1,4],sep="")
  }
  size_15_data = as.data.frame(rbindlist(serial_id_15,fill=TRUE))
  
  for (i in 1:length(rownames(size_15_data)))
  {
    if (size_15_data[i,3]>=32 & size_15_data[i,3]<=37)
    {
      loop[[i]] = "TL"
    }
    else if (size_15_data[i,3]>=54 & size_15_data[i,3]<=57)
    {
      loop[[i]] = "SL1"
    }
    else if (size_15_data[i,3]>=72 & size_15_data[i,3]<=77)
    {
      loop[[i]] = "SL2"
    }
    else if (size_15_data[i,3]>=88 & size_15_data[i,3]<=90)
    {
      loop[[i]] = "SL3"
    }
    else if (size_15_data[i,3]>=21 & size_15_data[i,3]<=31)
    {
      loop[[i]] = "R"
    }
    else if (size_15_data[i,3]>=38 & size_15_data[i,3]<=53)
    {
      loop[[i]] = "AR"
    }
    else if (size_15_data[i,3]>=63 & size_15_data[i,3]<=67)
    {
      loop[[i]] = "LR"
    }
    else
    {
      loop[[i]] = "NS"
    }
  }
  size_15_data["Structure"] = loop  
}
size_16 = subset(data_connection,data_connection[,ncol(data_connection)]==16)
if (length(row.names(size_16)) > 0)
{
  serial_id_list = unique(size_16[,4])
  serial_id_16 = NULL
  loop = NULL
  for (i in 1:length(serial_id_list))
  {
    serial_id_16[[i]] = subset(size_16,size_16[,4]==serial_id_list[[i]])
    serial_id_16[[i]]["Unique_ID"] = paste0(serial_id_16[[i]][1,2],"_",serial_id_16[[i]][1,3],"_",serial_id_16[[i]][1,4],sep="")
  }
  size_16_data = as.data.frame(rbindlist(serial_id_16,fill=TRUE))
  
  for (i in 1:length(rownames(size_16_data)))
  {
    if (size_16_data[i,3]>=32 & size_16_data[i,3]<=37)
    {
      loop[[i]] = "TL"
    }
    else if (size_16_data[i,3]>=54 & size_16_data[i,3]<=57)
    {
      loop[[i]] = "SL1"
    }
    else if (size_16_data[i,3]>=72 & size_16_data[i,3]<=77)
    {
      loop[[i]] = "SL2"
    }
    else if (size_16_data[i,3]>=88 & size_16_data[i,3]<=90)
    {
      loop[[i]] = "SL3"
    }
    else if (size_16_data[i,3]>=21 & size_16_data[i,3]<=31)
    {
      loop[[i]] = "R"
    }
    else if (size_16_data[i,3]>=38 & size_16_data[i,3]<=53)
    {
      loop[[i]] = "AR"
    }
    else if (size_16_data[i,3]>=63 & size_16_data[i,3]<=67)
    {
      loop[[i]] = "LR"
    }
    else
    {
      loop[[i]] = "NS"
    }
  }
  size_16_data["Structure"] = loop  
}
size_17 = subset(data_connection,data_connection[,ncol(data_connection)]==17)
if (length(row.names(size_17)) > 0)
{
  serial_id_list = unique(size_17[,4])
  serial_id_17 = NULL
  loop = NULL
  for (i in 1:length(serial_id_list))
  {
    serial_id_17[[i]] = subset(size_17,size_17[,4]==serial_id_list[[i]])
    serial_id_17[[i]]["Unique_ID"] = paste0(serial_id_17[[i]][1,2],"_",serial_id_17[[i]][1,3],"_",serial_id_17[[i]][1,4],sep="")
  }
  size_17_data = as.data.frame(rbindlist(serial_id_17,fill=TRUE))
  
  for (i in 1:length(rownames(size_17_data)))
  {
    if (size_17_data[i,3]>=32 & size_17_data[i,3]<=37)
    {
      loop[[i]] = "TL"
    }
    else if (size_17_data[i,3]>=54 & size_17_data[i,3]<=57)
    {
      loop[[i]] = "SL1"
    }
    else if (size_17_data[i,3]>=72 & size_17_data[i,3]<=77)
    {
      loop[[i]] = "SL2"
    }
    else if (size_17_data[i,3]>=88 & size_17_data[i,3]<=90)
    {
      loop[[i]] = "SL3"
    }
    else if (size_17_data[i,3]>=21 & size_17_data[i,3]<=31)
    {
      loop[[i]] = "R"
    }
    else if (size_17_data[i,3]>=38 & size_17_data[i,3]<=53)
    {
      loop[[i]] = "AR"
    }
    else if (size_17_data[i,3]>=63 & size_17_data[i,3]<=67)
    {
      loop[[i]] = "LR"
    }
    else
    {
      loop[[i]] = "NS"
    }
  }
  size_17_data["Structure"] = loop  
}

size_18 = subset(data_connection,data_connection[,ncol(data_connection)]==18)
if (length(row.names(size_18)) > 0)
{
  serial_id_list = unique(size_18[,4])
  serial_id_18 = NULL
  loop = NULL
  for (i in 1:length(serial_id_list))
  {
    serial_id_18[[i]] = subset(size_18,size_18[,4]==serial_id_list[[i]])
    serial_id_18[[i]]["Unique_ID"] = paste0(serial_id_18[[i]][1,2],"_",serial_id_18[[i]][1,3],"_",serial_id_18[[i]][1,4],sep="")
  }
  size_18_data = as.data.frame(rbindlist(serial_id_18,fill=TRUE))
  
  for (i in 1:length(rownames(size_18_data)))
  {
    if (size_18_data[i,3]>=32 & size_18_data[i,3]<=37)
    {
      loop[[i]] = "TL"
    }
    else if (size_18_data[i,3]>=54 & size_18_data[i,3]<=57)
    {
      loop[[i]] = "SL1"
    }
    else if (size_18_data[i,3]>=72 & size_18_data[i,3]<=77)
    {
      loop[[i]] = "SL2"
    }
    else if (size_18_data[i,3]>=88 & size_18_data[i,3]<=90)
    {
      loop[[i]] = "SL3"
    }
    else if (size_18_data[i,3]>=21 & size_18_data[i,3]<=31)
    {
      loop[[i]] = "R"
    }
    else if (size_18_data[i,3]>=38 & size_18_data[i,3]<=53)
    {
      loop[[i]] = "AR"
    }
    else if (size_18_data[i,3]>=63 & size_18_data[i,3]<=67)
    {
      loop[[i]] = "LR"
    }
    else
    {
      loop[[i]] = "NS"
    }
  }
  size_18_data["Structure"] = loop  
}
size_19 = subset(data_connection,data_connection[,ncol(data_connection)]==19)
if (length(row.names(size_19)) > 0)
{
  serial_id_list = unique(size_19[,4])
  serial_id_19 = NULL
  loop = NULL
  for (i in 1:length(serial_id_list))
  {
    serial_id_19[[i]] = subset(size_19,size_19[,4]==serial_id_list[[i]])
    serial_id_19[[i]]["Unique_ID"] = paste0(serial_id_19[[i]][1,2],"_",serial_id_19[[i]][1,3],"_",serial_id_19[[i]][1,4],sep="")
  }
  size_19_data = as.data.frame(rbindlist(serial_id_19,fill=TRUE))
  
  for (i in 1:length(rownames(size_19_data)))
  {
    if (size_19_data[i,3]>=32 & size_19_data[i,3]<=37)
    {
      loop[[i]] = "TL"
    }
    else if (size_19_data[i,3]>=54 & size_19_data[i,3]<=57)
    {
      loop[[i]] = "SL1"
    }
    else if (size_19_data[i,3]>=72 & size_19_data[i,3]<=77)
    {
      loop[[i]] = "SL2"
    }
    else if (size_19_data[i,3]>=88 & size_19_data[i,3]<=90)
    {
      loop[[i]] = "SL3"
    }
    else if (size_19_data[i,3]>=21 & size_19_data[i,3]<=31)
    {
      loop[[i]] = "R"
    }
    else if (size_19_data[i,3]>=38 & size_19_data[i,3]<=53)
    {
      loop[[i]] = "AR"
    }
    else if (size_19_data[i,3]>=63 & size_19_data[i,3]<=67)
    {
      loop[[i]] = "LR"
    }
    else
    {
      loop[[i]] = "NS"
    }
  }
  size_19_data["Structure"] = loop  
}
size_20 = subset(data_connection,data_connection[,ncol(data_connection)]==20)
if (length(row.names(size_20)) > 0)
{
  serial_id_list = unique(size_20[,4])
  serial_id_20 = NULL
  loop = NULL
  for (i in 1:length(serial_id_list))
  {
    serial_id_20[[i]] = subset(size_20,size_20[,4]==serial_id_list[[i]])
    serial_id_20[[i]]["Unique_ID"] = paste0(serial_id_20[[i]][1,2],"_",serial_id_20[[i]][1,3],"_",serial_id_20[[i]][1,4],sep="")
  }
  size_20_data = as.data.frame(rbindlist(serial_id_20,fill=TRUE))
  
  for (i in 1:length(rownames(size_20_data)))
  {
    if (size_20_data[i,3]>=32 & size_20_data[i,3]<=37)
    {
      loop[[i]] = "TL"
    }
    else if (size_20_data[i,3]>=54 & size_20_data[i,3]<=57)
    {
      loop[[i]] = "SL1"
    }
    else if (size_20_data[i,3]>=72 & size_20_data[i,3]<=77)
    {
      loop[[i]] = "SL2"
    }
    else if (size_20_data[i,3]>=88 & size_20_data[i,3]<=90)
    {
      loop[[i]] = "SL3"
    }
    else if (size_20_data[i,3]>=21 & size_20_data[i,3]<=31)
    {
      loop[[i]] = "R"
    }
    else if (size_20_data[i,3]>=38 & size_20_data[i,3]<=53)
    {
      loop[[i]] = "AR"
    }
    else if (size_20_data[i,3]>=63 & size_20_data[i,3]<=67)
    {
      loop[[i]] = "LR"
    }
    else
    {
      loop[[i]] = "NS"
    }
  }
  size_20_data["Structure"] = loop  
}

Serial_connection_data = rbind(size_20_data,size_19_data,size_18_data,size_17_data,size_15_data,size_14_data,size_13_data,size_12_data,size_11_data,size_10_data,size_9_data,size_8_data,size_7_data,size_5_data,size_4_data,size_3_data,size_2_data,size_1_data)

write.csv(Serial_connection_data,"Serial_connection_spacer_scaffold.csv",row.names = FALSE)



