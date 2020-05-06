#include <graphics.h>
#include <time.h>

int main()
{
    int i,j,p,q,gd,gm,G[700][700],t,s=16;
    srand(time(NULL));
	detectgraph(&gd,&gm);
	initgraph(&gd,&gm,NULL);
	setfillstyle(1,1);
	for(i=s;i<=1350;i=i+s)
    {
        for(j=s;j<=700;j=j+s)
        {
            if(rand()%20==1)
            {
                bar(i-s/2,j-s/2,i+s/2,j+s/2);
                G[i/s][j/s]=1;
            }
            else
            	G[i/s][j/s]=0;
        }
    }
	delay(2000);
	while(1)
	{
		for(i=s;i<=1350;i=i+s)
		{
			for(j=s;j<=700;j=j+s)
			{
				t=0;
				for(p=-1;p<=1;p++)
				{
					for(q=-1;q<=1;q++)
					{
						if(p==0&&q==0)
							continue;		
						if(getpixel(s+(i-s+(p*s)+s*(1350/s))%(s*(1350/s)),s+(j-s+(q*s)+s*(700/s))%(s*(700/s)))==1)
							t=t+1;	
					}
				}
				if(getpixel(i,j)==0)
				{
					if(t==2)
						G[i/s][j/s]=1;	
				}
				else
					G[i/s][j/s]=(getpixel(i,j)+1)%3;	
			}
		}
		for(i=s;i<=1350;i=i+s)
		{
			for(j=s;j<=700;j=j+s)
			{
				setfillstyle(1,G[i/s][j/s]);
				bar(i-s/2,j-s/2,i+s/2,j+s/2);
			}
		}
	}
	getch();
    closegraph();
    return 0;
}
