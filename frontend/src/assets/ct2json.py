ct=open('./ct.txt','r')
ct_cen=ct.readlines()
dataset={}
for line in ct_cen:
    li_list=line.rstrip().split(',')
    ds=li_list[0]
    ct=li_list[1]
    
    if ds not in dataset.keys():
        dataset[ds]=[]
        dataset[ds].append(ct)
    else:
        dataset[ds].append(ct)
        
print(dataset)