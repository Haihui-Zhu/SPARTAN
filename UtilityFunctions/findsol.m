function out=findsol(p,oof)

for i=1:length(oof)
    ptemp=p;
    ptemp(4,1)=p(4,1)-oof(i);
    outif(i,:)=roots(ptemp);
end

x=[0:1:30];
y=p(1,1).*x.^3+p(2,1).*x.^2+p(3,1).*x+p(4,1);
plot(x,y)
for i=1:length(outif(:,1))
    for j=1:length(outif(1,:))
        if isreal(outif(i,j))
            out(i,j)=outif(i,j);
        else
            out(i,j)=0;
        end
    end
end

end