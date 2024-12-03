clear all

session = '/media/data-313/Rat399-20190402/HPC/Rat399-20190402.xml';
SetCurrentSession(session);

dHPC = [25, 26, 30, 39];
vHPC = [22, 29, 31, 34, 35, 37, 38, 42];

s=GetSpikeTimes('all', 'output', 'full')
s(s(:,3)==0|s(:,3)==1,:)=[];
ClusNumb=unique(s(:,3));
Mf=[]

for i = dHPC
    g=s(find(s(:,2)==i),:)
    ClusNumb=unique(g(:,3));
    for j = ClusNumb
        t=s(find(s(:,3)==j),1)
        tBin=Bin(t.', 100)
        Mneuron=Accumulate(tBin)
        Mf=[Mf;Mneuron]
    end
end
    

    


