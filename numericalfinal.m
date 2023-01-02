f=uifigure('position',[20,30,800,800]);

listbox=uilistbox(f,'Items',{'Gauss elimination','Gauss seidle','Gauss jordan'},'position',[50 200 150 100]);
gobutton=uibutton(f,'text','GO','ButtonPushedFcn',@(gobutton,event) getmethod(gobutton,listbox,f));
global error_soln;
error_soln=[];
global flag;
flag=1;
global maxiterations;
maxiterations=12;
global counter;
counter=0;
function getmethod(gobutton,listbox,f)
listbox.Value
    if strcmp(listbox.Value,'Gauss elimination')==1
        delete(f.Children);
        label1=uilabel(f,'Text','Gauss Elimination','FontSize',20,'Position',[175 760 300 50]);
        label1=uilabel(f,'Text','Number of unknowns:','FontSize',12,'Position',[20 730 500 60]);
        textN=uieditfield(f,'FontSize',14,'Position',[150 750 50 15]);
        buttonN=uibutton(f,'FontSize',14,'text','ok','Position',[220 750 50 20],'ButtonPushedFcn',@(buttonN,event) getN(buttonN,textN,f));
    elseif strcmp(listbox.Value,'Gauss seidle')==1
        delete(f.Children);
        label1=uilabel(f,'Text','Gauss seilde','FontSize',20,'Position',[175 760 300 50]);
        label1=uilabel(f,'Text','Number of unknowns:','FontSize',12,'Position',[20 730 500 60]);
        textN=uieditfield(f,'FontSize',14,'Position',[150 750 50 15]);
        buttonN=uibutton(f,'FontSize',14,'text','ok','Position',[220 750 50 20],'ButtonPushedFcn',@(buttonN,event) getNsiedle(buttonN,textN,f));
    elseif strcmp(listbox.Value,'Gauss jordan')==1
        delete(f.Children);
        label1=uilabel(f,'Text','Gauss jordan','FontSize',20,'Position',[175 760 300 50]);
        label1=uilabel(f,'Text','Number of unknowns:','FontSize',12,'Position',[20 730 500 60]);
        textN=uieditfield(f,'FontSize',14,'Position',[150 750 50 15]);
        buttonN=uibutton(f,'FontSize',14,'text','ok','Position',[220 750 50 20],'ButtonPushedFcn',@(buttonN,event) getNjordan(buttonN,textN,f));
    end
end
function getNjordan(buttonN,textN,f)
    n=str2double(textN.Value);
    syms unknowns [n 1];
    stringat="";
    for(i=1:n)
        S=char(96+i);
        s=S;
        syms s;
        unknowns(i)=S;
        if(i~=n)
            stringat=strcat(stringat,S,",");
        else
            stringat=strcat(stringat,S);
        end
    end
    for i=1:n
        label2=uilabel(f,'Text',strcat('eqn:',num2str(i)),'FontSize',12,'Position',[20 730-i*20 145 60]);
        eqn(i)=uieditfield(f,'FontSize',14,'Position',[150 750-i*20 400 15])
    end
        calculate_btn=uibutton(f,'FontSize',14,'text','Calculate','Position',[220 750-(n+1)*30 100 20],'ButtonPushedFcn',@(calculate_btn,event) geteqnJordan(calculate_btn,f,n,eqn,stringat,unknowns));
end
function geteqnJordan(calculate_btn,f,n,eqn,stringat,unknowns)
    for(i=1:n)
        if i==1
            equations=[string(eqn(i).Value)];
        else
            equations=[equations,string(eqn(i).Value)];
        end
    end
    solution=[];
    matrixOfCoeff=zeros(n,n+1);
    for j=1:n
        str=strcat("@(",stringat,")",equations(j));
        fc=str2func(str);
        [b t]=coeffs(fc, unknowns);
        for i=1:n
          findIfZero=find(t==unknowns(i));
          if(isempty(findIfZero))
                matrixOfCoeff(j,i) = 0;
          else
            matrixOfCoeff(j,i) = b(findIfZero(1));
          end
        end
        findIfZero=find(t==1);
        if(isempty(findIfZero))
                matrixOfCoeff(j,n+1) = 0;
        else
            matrixOfCoeff(j,n+1) = b(findIfZero(1));
        end
    end
    timestart=datevec(datestr(now,'dd-mm-yyyy HH:MM:SS.FFF'));
    for i=1:n
        if matrixOfCoeff(i,i)==0
            fprintf("can't solve");
            break;
        else
            for j=n:-1:i+1
                if matrixOfCoeff(j,i)==0
                    continue;
                else
                    ratio=matrixOfCoeff(j,i)/matrixOfCoeff(i,i);
                    z=matrixOfCoeff(j,:)-ratio.*matrixOfCoeff(i,:);
                    matrixOfCoeff(j,:)=z;
                end
            end
            for v=i-1:-1:1
                if matrixOfCoeff(v,i)==0
                    continue;
                else
                    ratio=matrixOfCoeff(v,i)/matrixOfCoeff(i,i);
                    z=matrixOfCoeff(v,:)-ratio.*matrixOfCoeff(i,:);
                    matrixOfCoeff(v,:)=z;
                end
            end
        end
    end
    for i=1:n
       matrixOfCoeff(i,:)=matrixOfCoeff(i,:)./matrixOfCoeff(i,i);
       solution(i)=matrixOfCoeff(i,n+1);
    end
    timeend=datevec(datestr(now,'dd-mm-yyyy HH:MM:SS.FFF'));
    time=abs(etime(timestart,timeend))
    for(i=1:n)
        label=uilabel(f,'Text',strcat(char(96+i),':',num2str(solution(i))),'FontSize',12,'Position',[20 730-(n+3+i)*20 145 60]);
    end
    label=uilabel(f,'Text',strcat('time taken to calculate:',num2str(time)),'FontSize',12,'Position',[20 730-(2*n+3+1)*20 300 60]);
end
function getNsiedle(buttonN,textN,f)
    solution2=[];
    n=str2double(textN.Value);
    syms unknowns [n 1];
    stringat="";
    for(i=1:n)
        S=char(96+i);
        s=S;
        syms s;
        unknowns(i)=S;
        if(i~=n)
            stringat=strcat(stringat,S,",");
        else
            stringat=strcat(stringat,S);
        end
    end
    for i=1:n
        label2=uilabel(f,'Text',strcat('eqn:',num2str(i)),'FontSize',12,'Position',[20 730-i*20 145 60]);
        eqn(i)=uieditfield(f,'FontSize',14,'Position',[150 750-i*20 400 15])
    end
    for i=1:n
        label2=uilabel(f,'Text',strcat('init guess:',num2str(i)),'FontSize',12,'Position',[20 730-(n+i)*20 145 60]);
        initguess(i)=uieditfield(f,'FontSize',14,'Position',[150 750-(n+i)*20 400 15]);
    end
    label3=uilabel(f,'Text','max it.','FontSize',12,'Position',[20 730-(2*n+1)*20 145 60]);
    maxiterations=uieditfield(f,'FontSize',14,'Position',[150 750-(2*n+1)*20 400 15]);
    label3=uilabel(f,'Text','approx.','FontSize',12,'Position',[20 730-(2*n+2)*20 145 60]);
    app_error=uieditfield(f,'FontSize',14,'Position',[150 750-(2*n+2)*20 400 15]);
    calculate_btn=uibutton(f,'FontSize',14,'text','Calculate','Position',[220 750-(2*n+3)*22 100 20],'ButtonPushedFcn',@(calculate_btn,event) geteqnsielde(calculate_btn,f,n,eqn,stringat,unknowns,maxiterations,app_error,initguess));
end

function geteqnsielde(calculate_btn,f,n,eqn,stringat,unknowns,maxiterations,app_error,initguess)
    temp=""
    solution2=[];
    for i=1:n
        solution2(i)=str2num(initguess(i).Value);
    end
    maxiterations=str2num(maxiterations.Value);
    app_error=str2num(app_error.Value);
    flag=1;
    counter=0;
    for(i=1:n)
        if i==1
            equations=[string(eqn(i).Value)];
        else
            equations=[equations,string(eqn(i).Value)];
        end
    end
    matrixOfCoeff=zeros(n,n+1);
    for j=1:n
        str=strcat("@(",stringat,")",equations(j));
        fc=str2func(str);
        [b t]=coeffs(fc, unknowns);
        for i=1:n
          findIfZero=find(t==unknowns(i));
          if(isempty(findIfZero))
                matrixOfCoeff(j,i) = 0;
          else
            matrixOfCoeff(j,i) = b(findIfZero(1));
          end
        end
        findIfZero=find(t==1);
        if(isempty(findIfZero))
                matrixOfCoeff(j,n+1) = 0;
        else
            matrixOfCoeff(j,n+1) = b(findIfZero(1));
        end
    end
     timestart=datevec(datestr(now,'dd-mm-yyyy HH:MM:SS.FFF'));
    while(flag==1 && counter<maxiterations)
        counterTable=1;
        for i=1:n
            sum=0;
            for j=1:n
                if(j~=i)
                    sum=sum+matrixOfCoeff(i,j)*solution2(j);
                end
            end
            sum=sum+matrixOfCoeff(i,n+1);
            matrixOfCoeff(i,i);
            sum=-sum/matrixOfCoeff(i,i);
            err_soln(i)=abs((sum-solution2(i))/sum);
            solution2(i)=sum;
            tableElements(counter+1,counterTable)=solution2(i);
            tableElements(counter+1,counterTable+1)=err_soln(i);
            counterTable=counterTable+2
        end
        tableElements
        flag=0;
        for i=1:n
            if(err_soln(i)>app_error)
                flag=1;
            end
        end
        counter=counter+1;
        columnNames=[""]
    end
    j=1;
    for i=1:n
            columnNames(j)=string(unknowns(i));
            columnNames(j+1)=strcat("error:",string(unknowns(i)))
            j=j+2;
    end
    
    uit=uitable(f,'ColumnName',columnNames,'data',tableElements,'position',[20 750-(2*n+4)*20-300 750 300]);
    timeend=datevec(datestr(now,'dd-mm-yyyy HH:MM:SS.FFF'));
    time=abs(etime(timestart,timeend))
    for i=1:n
        label3=uilabel(f,'Text',strcat(char(96+i),':',num2str(solution2(i))),'FontSize',12,'Position',[20 750-(2*n+4+1)*20-300 1000 20]);
    end
    label4=uilabel(f,'Text',strcat('number of iterations',num2str(counter)),'FontSize',12,'Position',[20 750-(2*n+4+2)*20-300 1000 20]);
    label4=uilabel(f,'Text',strcat('time taken to calculate:',num2str(time)),'FontSize',12,'Position',[20 750-(2*n+4+3)*20-300 1000 20]);
    
end
function getN(buttonN,textN,f)
    n=str2double(textN.Value);
    syms unknowns [n 1];
    stringat="";
    for(i=1:n)
        S=char(96+i);
        s=S;
        syms s;
        unknowns(i)=S;
        if(i~=n)
            stringat=strcat(stringat,S,",");
        else
            stringat=strcat(stringat,S);
        end
    end
    for i=1:n
        label2=uilabel(f,'Text',strcat('eqn:',num2str(i)),'FontSize',12,'Position',[20 730-i*20 145 60]);
        eqn(i)=uieditfield(f,'FontSize',14,'Position',[150 750-i*20 400 15])
    end
        calculate_btn=uibutton(f,'FontSize',14,'text','Calculate','Position',[220 750-(n+1)*30 100 20],'ButtonPushedFcn',@(calculate_btn,event) geteqn(calculate_btn,f,n,eqn,stringat,unknowns));
end
function geteqn(calculate_btn,f,n,eqn,stringat,unknowns)
    for(i=1:n)
        if i==1
            equations=[string(eqn(i).Value)];
        else
            equations=[equations,string(eqn(i).Value)];
        end
    end
    solution=[];
    matrixOfCoeff=zeros(n,n+1);
    for j=1:n
        str=strcat("@(",stringat,")",equations(j));
        fc=str2func(str);
        [b t]=coeffs(fc, unknowns);
        for i=1:n
          findIfZero=find(t==unknowns(i));
          if(isempty(findIfZero))
                matrixOfCoeff(j,i) = 0;
          else
            matrixOfCoeff(j,i) = b(findIfZero(1));
          end
        end
        findIfZero=find(t==1);
        if(isempty(findIfZero))
                matrixOfCoeff(j,n+1) = 0;
        else
            matrixOfCoeff(j,n+1) = b(findIfZero(1));
        end
    end
    timestart=datevec(datestr(now,'dd-mm-yyyy HH:MM:SS.FFF'));
    for i=1:n-1
        if matrixOfCoeff(i,i)==0
            fprintf("can't solve");
            break;
        else
            for j=n:-1:i+1
                if matrixOfCoeff(j,i)==0
                    continue;
                else
                    ratio=matrixOfCoeff(j,i)/matrixOfCoeff(i,i);
                    z=matrixOfCoeff(j,:)-ratio.*matrixOfCoeff(i,:);
                    matrixOfCoeff(j,:)=z;
                end
            end
        end
    end
    solution(n)=-matrixOfCoeff(n,n+1)/matrixOfCoeff(n,n);
    for j=n-1:-1:1
        sum=0;
        for i=n:-1:j+1
            matrixOfCoeff(j,i)
            solution(i)
            sum=sum+matrixOfCoeff(j,i)*solution(i)
        end
        solution(j)=-(matrixOfCoeff(j,n+1)+sum)/matrixOfCoeff(j,j)
    end
    timeend=datevec(datestr(now,'dd-mm-yyyy HH:MM:SS.FFF'));
    time=abs(etime(timestart,timeend))
    for(i=1:n)
        label=uilabel(f,'Text',strcat(char(96+i),':',num2str(solution(i))),'FontSize',12,'Position',[20 730-(n+3+i)*20 145 60]);
    end
    label=uilabel(f,'Text',strcat('time taken to calculate:',num2str(time)),'FontSize',12,'Position',[20 730-(2*n+3+1)*20 300 60]);
end