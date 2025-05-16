int hindo()
{

    int sample[4][BUFSIZE]={0};
    int i=0;
    for(int i=0; i<BUFSIZE; i++)
    {
        
    for(int j=0; j<MAX_SEQ_NUM; j++)
    {
        if(g_motif[j][i]=='A')
        {
            sample[0][i]+=1;
        }
        if(g_motif[j][i]=='C')
        {
            sample[1][i]+=1;
        }    
        if(g_motif[j][i]=='G')
        {
            sample[2][i]+=1;
        }   
        if(g_motif[j][i]=='T')
        {
            sample[3][i]+=1;
        }
    }


    }
    for(int i=0; i<BUFSIZE; i++)
    {
        
    for(int j=0; j<MAX_SEQ_NUM; j++)
    {
        printf("%d\n",sample[j][i]);

    }

    }





    return 0;
}



/*

double p[4][BUFSIZE]={0};
    int A=0;

        for(int j=0; j<4; j++)
        {
            A=sample[j][0]++;        
        }

    for(int i=0; i<BUFSIZE; i++)
    {
        for(int j=0; j<4; j++)
        {
            p[j][i]=(1+sample[j][i])/(4+A); 
        }
    }

                            double q[4]={0};
                            int B=0;
                            for(int j=0; j<4; j++)
                            {
                                for(int i=0; i<BUFSIZE; i++)
                                {
                                    B=sample[j][i]++;
                                }
                            }
                            for(int j=0; j<4; j++)
                            {
                                for(int i=0; i<BUFSIZE; i++)
                                {
                                    int C=0;
                                    C=sample[j][i]++;
                                    q[j]=C/B;
                                }    
                            }

    double s[4][BUFSIZE]={0};
    for(int i=0; i<BUFSIZE; i++)
    {
        for(int j=0; j<4; j++)
        {
            s[j][i]=log((p[j][i])/(q[j])); 
        }
    }

    printf("%f\n",s[4][BUFSIZE]);

}



*/


