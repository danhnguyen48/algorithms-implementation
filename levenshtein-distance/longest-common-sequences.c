#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define size_c 200

#define EQUAL 1
#define UP 2
#define LEFT 3

int max(int m, int n, int i, int j);
int print(int i, int j);

int Lcs_len(char *str1, char *str2, int **char1, int **char2)
{
    unsigned long m=strlen(str1);
    unsigned long n=strlen(str2);
    int i,j;
    for(i=1; i<=m; i++){
        char1[i][0]=0;
    }
    for(j=1;j<=n;j++){
        char1[0][j]=0;
    }
    for(i=1;i<=m;i++){
        for(j=1;j<=n;j++){
            if(str1[i-1]==str2[j-1]){
                char1[i][j]=char1[i-1][j-1]+1;
                char2[i][j]=EQUAL;
            }
            else if (char1[i-1][j]>=char1[i][j-1]){
                char1[i][j]=char1[i-1][j];
                char2[i][j]=UP;
            }
            else{
                char1[i][j]=char1[i][j-1];
                char2[i][j]=LEFT;
            }
        }
    }
    return char1[m][n];
}

void Print_Lcs(char *str, int **b, long i, long j)
{
    if(i==0 || j==0)
        return;
    if(b[i][j]==EQUAL){
        Print_Lcs(str, b, i-1,j-1);
        printf("%c",str[i-1]);
    }
    else if (b[i][j]==UP)
        Print_Lcs(str, b, i-1, j);
    else
        Print_Lcs(str, b, i, j-1);
}

void Find_Lcs(char *str1,char *str2)
{
    int i,j, length;
    unsigned long len1=strlen(str1), len2=strlen(str2);
    int **c=(int**)malloc(sizeof(int*)*(len1+1));
    int **b=(int**)malloc(sizeof(int*)*(len1+1));
    for(i=0;i<=len1;i++){
        c[i]=(int*)malloc(sizeof(int)*(len2+1));
        b[i] = (int *)malloc(sizeof(int) * (len2 + 1));
    }
    for ( i = 0; i<= len1; i++)
    for( j = 0; j <= len2; j++){
        c[i][j] = 0;
        b[i][j] = 0;
    }
    length = Lcs_len(str1, str2, c, b);
    printf("The number of the Longest-Common-Subsequence is %d\n", length);
    printf("The Longest-Common-Subsequence is: ");
    Print_Lcs(str1, b, len1, len2);
    printf("\n");
    for ( i = 0; i <= len1; i++){
        free(c[i]);
        free(b[i]);
    }
    free(c);
    free(b);

}

int main(int argc, const char * argv[]) {
    char X[size_c],Y[size_c];
    printf("please enter your characters:");
    scanf("%s",X);
    
    while(strlen(X) > 200){
        printf("what you input is too long, please try again");
        scanf("%s\n",X);
    }
    printf("please enter your characters:");
    scanf("%s",Y);
    while(strlen(Y) > 200){
    printf("what you input is too long, please try again");
    scanf("%s",Y);
        }
    Find_Lcs(X,Y);
    return 0;
}
