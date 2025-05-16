#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define BUFSIZE 1024 //ファイルから読み込む一行の最大文字数
#define MAX_SEQ_NUM 30 //一つの転写因子に対して与えられる結合部位配列の最大数
#define MAX_GENE_NUM 8 /*与えられるプロモータ領域の最大遺伝子数*/
#define SUUJI 4

char g_motif[MAX_SEQ_NUM][BUFSIZE]; //転写因子の結合部位配列を保存する配列


struct promoter{
  char name[BUFSIZE];
  char seq[BUFSIZE];
}g_pro[MAX_GENE_NUM]; //遺伝子のプロモータ領域を保存する構造体


int sample[SUUJI][BUFSIZE] = {0};
double p[SUUJI][BUFSIZE] = {0};
double q[SUUJI] = {0};
double s[SUUJI][BUFSIZE] = {0};
double hit[BUFSIZE] = {0};

int seq_num, gene_num, motif_len;

//グローバル変数はローカル変数と区別するため、全部大文字にするかg_を先頭につけるのが一般的

int read_multi_seq(char* filename){
  seq_num = 0;
  char buffer[BUFSIZE];
  FILE *fp = fopen(filename,"r");

  if(fp == NULL){
    printf("motif_region_file open error.\n");
    exit(1); //ファイルが開けなかった場合プログラムを終了
  }

  while(fscanf(fp, "%s", buffer) != EOF){ //プログラムから一行ずつ読み込む
    if(buffer[strlen(buffer)-1]=='\n'){
      buffer[strlen(buffer)-1]='\0'; //改行を切り落とす
    }
    strcpy(g_motif[seq_num],buffer); //結合部位配列を保存
    seq_num++;
  }
  return seq_num;
}

int read_promoter(char *filename){
  int gene_num = 0;  
  char buffer[BUFSIZE];
  FILE *fp = fopen(filename,"r");

  if(fp == NULL){
    printf("scorefile open error.\n");
    exit(1);
  }

  while(fscanf(fp, "%s", buffer) != EOF){
    if(buffer[strlen(buffer)-1]=='\n'){
      buffer[strlen(buffer)-1]='\0';
    }
    
    if(buffer[0]=='>'){
      strcpy(g_pro[gene_num].name,buffer+1); 
    }else{
      strcpy(g_pro[gene_num].seq,buffer);
      gene_num++;
    }    
  }
  return gene_num;
}


    int c_sample()
    {
    int i=0;
    for(int i=0; i<motif_len; i++)
    {
        
    for(int j=0; j<seq_num; j++)
    {
        if(g_motif[j][i]=='A')
        {
            sample[0][i]+=1;
        }
        else if(g_motif[j][i]=='C')
        {
            sample[1][i]+=1;
        }    
        else if(g_motif[j][i]=='G')
        {
            sample[2][i]+=1;
        }   
        else if(g_motif[j][i]=='T')
        {
            sample[3][i]+=1;
        }

    }



    }

    printf("Position\tA\tC\tG\tT\n");
    for (int i = 0; i < motif_len; i++) 
    {
        printf("%d\t\t%d\t%d\t%d\t%d\n",
            i+1,
            sample[0][i],
            sample[1][i],
            sample[2][i],
            sample[3][i]
        );
    }

    }






    void c_p()
    {
    int A=seq_num+SUUJI;

    for(int i=0; i<motif_len; i++)
    {
        for(int j=0; j<SUUJI; j++)
        {
            p[j][i]=(1.0+sample[j][i])/A; 
        }
    }

    printf("\n");
    printf("Position\tA\t\tC\t\tG\t\tT\n");
    for (int i = 0; i < motif_len; i++) 
    {
        printf("%d\t\t%f\t%f\t%f\t%f\n",
            i+1,
            p[0][i], 
            p[1][i], 
            p[2][i], 
            p[3][i]  
        );
    }

    }


    void c_q()
    {
    double B=0;
    double b[SUUJI]={7519429, 4637676, 4637676, 7519429};

        for(int i=0; i<SUUJI; i++)
        {
            B+=b[i];
        }
    
    for(int j=0; j<SUUJI; j++)
    {
    
        q[j]=b[j]/B;
        
    }
    

    }


    void c_s()
    {
    for(int i=0; i<motif_len; i++)
    {
        for(int j=0; j<SUUJI; j++)
        {
            s[j][i]=log((p[j][i])/(q[j])); 
        }
    }

    printf("\n");
    printf("Position\tA\t\tC\t\tG\t\tT\n");
    for (int i = 0; i < motif_len; i++) 
    {
        printf("%d\t\t%f\t%f\t%f\t%f\n",
            i+1,
            s[0][i], 
            s[1][i], 
            s[2][i], 
            s[3][i]  
              );
    }
    
    }



    void c_hit()
    {
   
    
    for(int i=0; i<gene_num; i++)
    {
      int len= strlen(g_pro[i].seq);
      double MAX=-10000;
      int best=0;

      for(int j=0; j <= len - motif_len; j++)
      {
        double total=0;
                  for(int g=0; g<motif_len; g++)
                 {
                  int n=300;

                  if(g_pro[i].seq[j+g]=='A')
                  {
                    n=0;
                  }

                  else if(g_pro[i].seq[j+g]=='C')
                  {
                    n=1;
                  }

                  else if(g_pro[i].seq[j+g]=='G')
                  {
                    n=2;
                  }

                  else if(g_pro[i].seq[j+g]=='T')
                  {
                    n=3;
                  }

                  if(n!=300)
                  {
                  total+=s[n][g];
                  }                     
                 }
          
          if(total>MAX)
          {
            MAX=total;
            best=j+1;
          }

      }

      printf("\n");
      printf("Gene %s: Max Score = %f at position %d", g_pro[i].name, MAX, best);


    }
      printf("\n\n");
    }



int main(int argc, char* argv[]){
  seq_num = read_multi_seq(argv[1]); //１番目の引数で指定した転写因子の複数の結合部位配列を読み込む

//  printf("motif region:\n");
  for(int i = 0; i < seq_num; i++){
//    printf("%s\n",g_motif[i]); //読み込んだ転写因子の結合部位配列を表示
  }
  printf("\n");

  gene_num = read_promoter(argv[2]);  //２番目の引数で指定した遺伝子のプロモータ領域を読み込む
  
//  printf("promoter_sequence:\n");
  for(int i = 0; i < gene_num; i++){
//    printf(">%s\n", g_pro[i].name); //読み込んだプロモータ領域を表示
    //printf("%s\n", g_pro[i].seq);
  }

  motif_len=strlen(g_motif[0]);

    c_sample();
    c_p();
    c_q();
    c_s();

    c_hit();

  return 0;
}
    
