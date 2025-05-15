#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define BUFSIZE 1024 //ファイルから読み込む一行の最大文字数
#define MAX_SEQ_NUM 30 //一つの転写因子に対して与えられる結合部位配列の最大数
#define MAX_GENE_NUM 8 /*与えられるプロモータ領域の最大遺伝子数*/

char g_motif[MAX_SEQ_NUM][BUFSIZE]; //転写因子の結合部位配列を保存する配列


struct promoter{
  char name[BUFSIZE];
  char seq[BUFSIZE];
}g_pro[MAX_GENE_NUM]; //遺伝子のプロモータ領域を保存する構造体

//グローバル変数はローカル変数と区別するため、全部大文字にするかg_を先頭につけるのが一般的

int read_multi_seq(char* filename){
  int seq_num = 0;
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



int main(int argc, char* argv[]){
  int seq_num = read_multi_seq(argv[1]); //１番目の引数で指定した転写因子の複数の結合部位配列を読み込む

//  printf("motif region:\n");
  for(int i = 0; i < seq_num; i++){
//    printf("%s\n",g_motif[i]); //読み込んだ転写因子の結合部位配列を表示
  }
  printf("\n");

  int gene_num = read_promoter(argv[2]);  //２番目の引数で指定した遺伝子のプロモータ領域を読み込む
  
//  printf("promoter_sequence:\n");
  for(int i = 0; i < gene_num; i++){
//    printf(">%s\n", g_pro[i].name); //読み込んだプロモータ領域を表示
    //printf("%s\n", g_pro[i].seq);
  }

    int sample[4][BUFSIZE]={0};
    int motif_len = strlen(g_motif[0]);
    int i=0;
    for(int i=0; i<motif_len; i++)
    {
        
    for(int j=0; j<MAX_SEQ_NUM; j++)
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





    double p[4][BUFSIZE]={0};
    int A=0;

    A=seq_num+4;

    for(int i=0; i<motif_len; i++)
    {
        for(int j=0; j<4; j++)
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

    double q[4]={0};
    double B=0;
    double b[4]={7519429, 4637676, 4637676, 7519429};

        for(int i=0; i<4; i++)
        {
            B+=b[i];
        }
    
    for(int j=0; j<4; j++)
    {
    
        q[j]=b[j]/B;
        
    }

    double s[4][BUFSIZE]={0};
    for(int i=0; i<motif_len; i++)
    {
        for(int j=0; j<4; j++)
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


    /*
    double hit[BUFSIZE];


    for(int i=0; i<motif_len; i++)
    {
        if(g_pro[i]=='A')
        {
            hit[]+=s[0][i]
        }

        else if(=='C')
        {
            hit[]+=s[1][i]
        }

        else if(=='G')
        {
            hit[]+=s[2][i]
        }

        else if(=='T')
        {
            hit[]+=s[3][i]
        }


    }


*/
    
  return 0;
}



    
