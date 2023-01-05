#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#define SIZE 4096

struct DataItem {
   int data;   
   int key;
};

struct DataItem* hashArray[SIZE]; 
struct DataItem* dummyItem;
struct DataItem* item;

int hashCode(int key) {
   return key % SIZE;
}

struct DataItem *search(int key) {
   //get the hash 
   int hashIndex = hashCode(key);  
	
   //move in array until an empty 
   while(hashArray[hashIndex] != NULL) {
	
      if(hashArray[hashIndex]->key == key)
         return hashArray[hashIndex]; 
			
      //go to next cell
      ++hashIndex;
		
      //wrap around the table
      hashIndex %= SIZE;
   }        
	
   return NULL;        
}

void insert(int key,int data) {

   struct DataItem *item = (struct DataItem*) malloc(sizeof(struct DataItem));
   item->data = data;  
   item->key = key;

   //get the hash 
   int hashIndex = hashCode(key);

   //move in array until an empty or deleted cell
   while(hashArray[hashIndex] != NULL && hashArray[hashIndex]->key != -1) {
      //go to next cell
      ++hashIndex;
		
      //wrap around the table
      hashIndex %= SIZE;
   }
	
   hashArray[hashIndex] = item;
}

struct DataItem* delete(struct DataItem* item) {
   int key = item->key;

   //get the hash 
   int hashIndex = hashCode(key);

   //move in array until an empty
   while(hashArray[hashIndex] != NULL) {
	
      if(hashArray[hashIndex]->key == key) {
         struct DataItem* temp = hashArray[hashIndex]; 
			
         //assign a dummy item at deleted position
         hashArray[hashIndex] = dummyItem; 
         return temp;
      }
		
      //go to next cell
      ++hashIndex;
		
      //wrap around the table
      hashIndex %= SIZE;
   }      
	
   return NULL;        
}

void display() {
   int i = 0;
	
   for(i = 0; i<SIZE; i++) {
	
      if(hashArray[i] != NULL)
         printf(" (%d,%d)",hashArray[i]->key,hashArray[i]->data);
      else
         printf(" ~~ ");
   }
	
   printf("\n");
}

int kmer_value(char kmer[])
{
    return pow(4,5)*(u_int8_t)(kmer[0]-'0')+pow(4,4)*(u_int8_t)(kmer[1]-'0')+pow(4,3)*(u_int8_t)(kmer[2]-'0')+pow(4,2)*(u_int8_t)(kmer[3]-'0')+pow(4,1)*(u_int8_t)(kmer[4]-'0')+pow(4,0)*(u_int8_t)(kmer[5]-'0');
}

void mm_kmermap_seq(int qlen, u_int8_t *input_seq[2], int output_qlen, u_int8_t *output_seq[2]) {
   //dummyItem = (struct DataItem*) malloc(sizeof(struct DataItem));
   //dummyItem->data = -1;  
   //dummyItem->key = -1;

   //printf("kmer value: %d\n", kmer_value("001000"));
   int qlevel = 3;
   FILE* ptr;
   char aline[20];
   ptr = fopen("quantizedmap.txt","r");
   int count = 0;
   printf("\n\nStarting the while loop: ptr = %p; aline = %s", ptr, aline);
   while(fscanf(ptr,"%s\n", aline)==1){
      //printf("%s\n",aline);
      if(count >= 4){
         int akey;
         int val;
         //printf("%s\t",aline);
         if (count%4 == 0){
            akey = kmer_value(aline);
         } else if ((count-qlevel+1)%4 == 0){
            val = *aline - '0';
            //printf("%d ", val);
         }
         if ((count+1)%4 == 0){
            insert(akey,val);
         }
      }
      count++;
   }
   fclose(ptr);
   
   //u_int8_t *output_seq[2];
   //int kmer = 6;
   //int output_qlen = qlen - kmer + 1;
   //output_seq[0] = (u_int8_t*)malloc(output_qlen * 2);
	//output_seq[1] = output_seq[0] + output_qlen;
   for (int i = 0; i < output_qlen; ++i){
      char akey0[6], akey1[6];
      for (int j = 0; j < 6; ++j){
         akey0[j] = input_seq[0][i+j] + '0';
         akey1[j] = input_seq[1][i+j] + '0';
      }
      int item0 = search(kmer_value(akey0))->data;
      int item1 = search(kmer_value(akey1))->data;
      printf("\ni=%d, %d, %d", i, item0, item1);
      output_seq[0][i] = item0;
      output_seq[1][i] = item1;
   }
   /*
   for (int i = 0; i < output_qlen; ++i){
      if (i == 0){
         printf("\noutput_seq0: ");
      }
      printf("%d", output_seq[0][i]);
   }
   for (int i = 0; i < output_qlen; ++i){
      if (i == 0){
         printf("\noutput_seq1: ");
      }
      printf("%d", output_seq[1][i]);
   }*/
return;
   /*
   insert(kmer_value("001001"), 0);
   insert(kmer_value("001010"), 1);
   insert(kmer_value("001110"), 2);
   insert(kmer_value("001010"), 1);
   insert(kmer_value("001002"), 2);
   insert(kmer_value("001003"), 2);
   insert(kmer_value("001200"), 0);
   insert(kmer_value("001231"), 1);
   insert(kmer_value("001111"), 0);
   
   display();
   item = search(kmer_value("001000"));

   if(item != NULL) {
      printf("Element found: %d\n", item->data);
   } else {
      printf("Element not found\n");
   }

   //delete(item);
   item = search(kmer_value("001302"));

   if(item != NULL) {
      printf("Element found: %d\n", item->data);
   } else {
      printf("Element not found\n");
   }*/
}


int main(){
   u_int8_t *seq[2];
   int qlen = 10;
   seq[0] = (u_int8_t*)malloc(qlen * 2);
   seq[1] = seq[0] + qlen;
   for (int i = 0; i < qlen; ++i){
      seq[0][i] = i % 4;
      seq[1][qlen - i - 1] = 3 - seq[0][i];
   }

   for(int i=0; i<qlen; ++i){
      if (i == 0){
         printf("\nseq0: ");
      }
      printf("%d", seq[0][i]);
   }

   for(int i=0; i<qlen; ++i){
      if (i == 0){
         printf("\nseq1: ");
      }
      printf("%d", seq[1][i]);
   }

   int out_qlen = qlen-5;
   u_int8_t *out_seq[2];
   out_seq[0] = (u_int8_t*)malloc(out_qlen * 2);
   out_seq[1] = out_seq[0] + out_qlen;
   
   mm_kmermap_seq(qlen,seq, out_qlen, out_seq);

   for(int i=0; i<out_qlen; ++i){
      if (i == 0){
         printf("\nout_seq0: ");
      }
      printf("%d", out_seq[0][i]);
   }

   for(int i=0; i<out_qlen; ++i){
      if (i == 0){
         printf("\nout_seq1: ");
      }
      printf("%d", out_seq[1][i]);
   }
   
   return 0;
}