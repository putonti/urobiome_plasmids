# written as part of work entitled "Plasmids of the Urinary Microbiota"


import os

# PARAMETERS TO SET BEFORE RUNNING
path_to_mash='mypath' ###CHANGE### path where mash is. If in path, remove this variable and future reference to it
list_of_files='my_list.txt' ###CHANGE### This is a list of the files to be analyzed
mash_threshold=0.01 ###CHANGE IF YOU WANT TO ALLOW MORE/LESS SEQUENCE SIMILARITY###

out_path='/mypath/' ###CHANGE### where output of mash and cluster sequences will be written
out_file='mash_results.txt' ###CHANGE### name of file where you want output


# RUN MASH CALCULATIONS

#get file names
with open(list_of_files) as f:
    d=f.readlines()
data=[i.strip() for i in d]

# make mash calls
command=path_to_mash + 'mash dist '
c1=0
for i in data:
    c2 = 0
    for j in data[c1:]:
        call=command+i+' '+j+' | tee -a '+out_path+out_file
        os.system(call)
        c2+=1
    c1+=1


# PARSE MASH RESULTS

with open(list_of_files) as f:
    d=f.readlines()
data=[i.strip() for i in d]
n=len(data)

# initialize 2D array
mash_values=[]
for i in range(n):
    mash_values.append([])
    for j in range(n):
        mash_values[i].append(0)

with open(out_path+'mash_results.txt','r') as f:
    x=f.readlines()    
x=[i.strip() for i in x]

for i in x:
    v=i.split('\t')
    s1=v[0]
    s2=v[1]
    m_dist=float(v[2])
    p_value=float(v[3])
    n_match=int(v[4][:v[4].find('/')])

    x=data.index(s1)
    y=data.index(s2)
    mash_values[x][y]=m_dist
    if x!=y:
        mash_values[y][x]=m_dist

# write out matrix as a csv file
with open(out_path+'mash_matrix.csv','w') as f:
    f.write(','.join(data))
    f.write('\n')
    for i in range(n):
        for j in range(n):
            f.write(str(mash_values[i][j])+',')
        f.write('\n')



# FIND CLUSTERS
clusters=dict()
counter=0
for i in range(n):
    for j in range(i+1,n):
        if mash_values[i][j]<mash_threshold:
            # i matches j
            # already in a cluster?
            found=False
            for k in clusters:
                if i in clusters[k] or j in clusters[k]:
                    clusters[k].add(i)
                    clusters[k].add(j)
                    found=True
            if found==False:
                clusters[counter]=set()
                clusters[counter].add(i)
                clusters[counter].add(j)
                counter+=1

# write out clusters to file (in folder called clusters in out_path)
for i in clusters:
    with open(out_path+'clusters/cluster_'+str(i)+'.fasta','w') as outf:
        for j in clusters[i]:
            with open(data[j],'r') as f:
                x=f.readlines()
            outf.writelines(x)


