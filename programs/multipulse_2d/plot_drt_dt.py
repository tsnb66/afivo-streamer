import matplotlib.pyplot as plt



#List of files
fileList = ["./multipulse_"+x+"_output.txt" for x in ["simple", "complex"]]
fileListdt = []
keys = ["cfl", "diff", "drt", "chem"]

#Read each file and get a dict of all their dts
for i,file in enumerate(fileList):
    fileListdt.append({x: [] for x in keys})
    flist = [x.strip() for x in open(file,"r").readlines()]
    lines = [list(map(float,x.split()[1:5])) for x in flist if x.startswith('dt')]
    
    #Print the cfl times
    for j in range(len(keys)):
        plt.figure(j)
        plt.title(keys[j])
        plt.plot([x[j] for x in lines], label=file)


plt.legend()
plt.show()




