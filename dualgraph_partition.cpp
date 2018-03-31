#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <ctime>
#include <stack>
#include<string.h>
#include<stdlib.h>
using namespace std;
//node
struct edge
{ int v1;
  int v2;
};
struct node
{
    int vertex;
	int binary;
	float rel;
    node * next;
	
};
int nods =0;
int motif_counter= 0;
int numreg_blocks = 0;
int numpseudo_blocks = 0;
int edgs =0;
int dim = 0;
int iterat = 0;
int term = 0;
stack<edge> mystack;
stack<edge> mystack1;
bool * Visited;
int *d_cmpnt; 
node * glovalnode;
int  weight[1000][1000]; //weight[100][100] weight[200][200]
ifstream infile;
ofstream outfile;
 bool cycles=false;// determine if a graph has cycles.
 
 
class Graph {
private:
	int n;
    int e; // number of edges
	float A;
	int B;
	int b_num1;
	int b_low1;
	int * b_num;
	int * b_low;
	int D;
	edge x;
	int terminal;
	long counter;
	int * temp;
	
	bool * T;
	int * L;
	int * father;
	node * w;
    edge* edge_irrel;
    int * Ls;
    int * Lt;
public:  
	 node * headnodes;
	 node * headnodes_A;
	 Graph (int nods, int edgs, int dim, int term) 
	// construtor
	{  
		//infile.open("C:/users/petingi/desktop/RnaDv/V7adj.txt");
	    terminal=term;
        counter;
		n=nods;
		b_num1 = 0;
		b_low1 = 0;
		b_low = new int [n];
		b_num = new int [n];
		
		e = edgs;
		D= dim;
		counter=0;
        headnodes= new node [n];
		//int temp[200];			//temporary data structure which is needed to 
		T = new bool[n];
		L= new int[n];
		father = new int[n];
		d_cmpnt = new int [n];
        temp=new int[n];
        edge_irrel = new edge [e];
        Ls = new int [n];
        Lt = new int [n];
        //headnodes_A= new node [n];
		// cout << "ok until here 2" << endl;
		// headnodes is an array of nodes.
        for (int i=0; i < n; i++)
         {  
			 headnodes[i].vertex=i;
             headnodes[i].next=0;
			 b_num1 = 1;
		     b_low1 = 0;
		     b_low [i] = 0;
		     b_num [i] = 0;
			 d_cmpnt [i] = 0;
         }
		// cout << "ok until here 3" << endl;
	}
	  void init () 
	// construtor
	{  
		//infile.open("C:/users/petingi/desktop/RnaDv/V7adj.txt");
	    //terminal=term;
        counter=0;
		n=nods;
		b_num1 = 0;
		numreg_blocks = 0;
        numpseudo_blocks = 0;
		//b_low = new int [n];
		//b_num = new int [n];
		
		//e = edgs;
		//D= dim;
		//counter=0;
        //headnodes= new node [n];
		//int temp[200];			//temporary data structure which is needed to 
		//T = new bool[n];
		//L= new int[n];
		//father = new int[n];
		//d_cmpnt = new int [n];
        //temp=new int[n];
        //edge_irrel = new edge [e];
        //Ls = new int [n];
        //Lt = new int [n];
        //headnodes_A= new node [n];
		//cout << "ok until here 2: " << n << endl;
		// headnodes is an array of nodes.
        for (int i=0; i < n; i++)
         {  
			 headnodes[i].vertex=i;
             headnodes[i].next=0;
		     b_num1 = 1;
		     b_low1 = 0;
		     b_low [i] = 0;
		     b_num [i] = 0;
			 d_cmpnt [i] = 0;
         } 
		    while ( !mystack.empty()) 
			  mystack.pop();
			while ( !mystack1.empty()) 
			  mystack1.pop(); 
	}
	 ~Graph ()
	 {  delete [ ] headnodes;
		//int temp[200];			//temporary data structure which is needed to 
		delete [] T;
		delete [] L;
		delete [] father;
        delete [] edge_irrel;
        delete [] Ls;
        delete [] Lt;
	 }
	 void Del ( )
	 {  delete [] headnodes;
	   //delete [] headnodes_A;
		delete [] T;
		delete [] L;
		delete [] father;
        delete [] edge_irrel;
        delete [] Ls;
        delete [] Lt;
		//delete [][] weight;
	 }
	
void setdiam (int diameter)
{ D=diameter;
}
	float factoring (int flag)
	{   
		int value;
		++counter;
		
	    Dijkstra_B (value, x);
		//cout << endl << "value: " << value << " " << x.v1 << " " << x.v2;
		if (value == 1) return 1;
		else if (value==0) return 0;
		else { // value = 2
			
		node *temp1, *nextn;
	    temp1 = &headnodes[x.v1];//delete edge
		//cout << endl << "vertex 1 " << x.v1	<< "vertex 2" << x.v2;
		while(temp1->vertex!= x.v2)
			{	
				temp1 = temp1->next; // tranfer the address of 'temp->next' to 'temp'
			}   
			
			temp1->binary = 1;
			temp1->rel = 1;
			float second_call= factoring (0);
			
			temp1->binary = 0; //delete edge
			temp1->rel = 0.5;
			//cout << endl << "vertex: " << temp1->vertex;
			float first_call = factoring (1);
			temp1->binary = 1;
			temp1->rel = 0.5;
			
			return (0.5 * first_call + 0.5 * second_call);
		}
        
	} 
	float factoring_irrel (int flag)
	{   
		int value, y, r, irr;
		edge x1;
		++counter;
		int inter_counter=counter;
		int edge_counter = 0;
		edge* edge_irrel_1,*edge_irrel_2;
		//cout << endl << "counter: " << counter;
        node *temp1, *temp2;
		Dijkstra_B (value, x1);
		
		if (value == 1) { 
			
			  return (1);
		}
		else if (value==0) {
            
			return 0;
		}
		else {  if (flag)  { 
                edge_irrel_2=new edge [2*e];
			    irr=Dltirrlvnt(edge_irrel_2, edge_counter);
				if (!irr) delete [ ] edge_irrel_2;
				if (irr) 
				{   edge_irrel_1 = edge_irrel_2;
					Dijkstra_B (value, x1);
				  //cout << endl << "+++++++++++++++++++++++++++";
                  //cout << endl << "cuando hay irrel edge " << x1.v1 << "-" << x1.v2;
				  if (value < 2){
                       for (int j=0; j < edge_counter ; j++)
			             { y = edge_irrel_1[j].v1;
			               r = edge_irrel_1[j].v2;
                           temp1 = &headnodes[y];//delete edge
		                    while(temp1->vertex!= r)
			                {	
				              temp1 = temp1->next; // tranfer the address of 'temp->next' to 'temp'
			                }  
					         temp1->binary = 1;
					         temp1->rel = 0.5;
					         //cout << endl << "======= inderting irrel ======";
					         //cout << endl << "edge irrelevante: (" << y <<", "<< r <<")";
                             //cout << endl << "======= inderting irrel ======";
					   }// for*/
					   //display_val();
					   //cout << endl << "counter: " << counter;
					   delete [ ] edge_irrel_2;
                       if (value ==1){
					       //cout << endl << " value irr: " << 1;
			               //display_val ()
						   return (1);
					   }
					   else {//cout << endl << " value irr: " << 1;
			                 //display_val (); 
							 return (0);
							 }
				//}
				  }//value
				}
		}
		//cout << endl << "============== edge " << x.v1 << "-" << x.v2;
	    temp1 = &headnodes[x1.v1];//delete edge
		while(temp1->vertex!= x1.v2)
			{	
				temp1 = temp1->next; // tranfer the address of 'temp->next' to 'temp'
			}  
		
			temp1->binary = 1;
			temp1->rel = 1;
        temp2 = &headnodes[x1.v2];//delete edge
		while(temp2->vertex!= x1.v1)
			{	
				temp2 = temp2->next; // tranfer the address of 'temp->next' to 'temp'
			}  //delete edge
			//temp2->rel = 0.5;
			
		
			temp2->binary = 1;
			temp2->rel = 1;
			//cout << endl << "primero: " << temp2->vertex;*/
			//cout << endl << "segundo: " << temp1->vertex;
			float second_call= factoring_irrel (0);
			
			temp1->binary = 0; //delete edge
			temp1->rel = 0.5;
            temp2->binary = 0; //delete edge
			temp2->rel = 0.5;
            //temp2->binar
			//cout << endl << "vertex: " << temp1->vertex;
			float first_call = factoring_irrel (1);
			
              temp1->binary = 1; //delete edge
		      temp1->rel = 0.5;
              temp2->binary = 1; //delete edge
			  temp2->rel = 0.5;
			  if (flag && irr) // put back irrelevant edges
			  { //cout << endl << "============== edge " << x1.v1 << "-" << x1.v2;
				  
				  for (int j=0; j < edge_counter ; j++)
			        { y = edge_irrel_1[j].v1;
			          r = edge_irrel_1[j].v2;
                      temp1 = &headnodes[y];//delete edge
		              while(temp1->vertex!= r)
			           {	
				         temp1 = temp1->next; // tranfer the address of 'temp->next' to 'temp'
			           }  
					  temp1->binary = 1;
					  temp1->rel = 0.5;
				  }
				  //delete [ ] edge_irrel_2;
			  }
			return (0.5 * first_call + 0.5 * second_call);
		}
        
	} 
float factoring_irrel_C (int flag)
	{   
		int value, y, r, irr, u1;
		++counter;
		int edge_counter = 0;
		edge* edge_irrel_1,*edge_irrel_2;;
		edge x1;
        node *temp1, *temp2;
		Dijkstra_B (value, x1);
		if (value == 1) return (1);
		else if (value==0) return (0);
		else {  if (flag)  {  
			    edge_irrel_2 = new edge [2*e];
	            irr=Dltirrlvnt_C(edge_irrel_2,edge_counter);
				if (!irr) delete [] edge_irrel_2;
                if (irr){
				edge_irrel_1= edge_irrel_2;
                Dijkstra_B (value, x1);
				if (value < 2) {
					//cout << endl << "retorno counter: " << edge_irrel;
                              for (int j=0; j < edge_counter ; j++)
			                  { y = edge_irrel_1[j].v1;
			                    r = edge_irrel_1[j].v2;
                                temp1 = &headnodes[y];//delete edge
		                        while(temp1->vertex!= r)
			                  {	
				                temp1 = temp1->next; // tranfer the address of 'temp->next' to 'temp'
			                  }  
					            temp1->binary = 1;
					            temp1->rel = 0.5;
							  }
							 delete [ ] edge_irrel_2;
                             return (value);
				  }
					
				}
		}
		//cout << endl << "============== edge " << x.v1 << "-" << x.v2;
	    temp1 = &headnodes[x1.v1];//delete edge
		temp2 = &headnodes[x1.v2];
		while(temp1->vertex!= x1.v2)
			{	
				temp1 = temp1->next; // tranfer the address of 'temp->next' to 'temp'
			}  
        while(temp2->vertex!= x1.v1)
			{	
				temp2 = temp2->next; // tranfer the address of 'temp->next' to 'temp'
			}  //delete edge
			//temp2->rel = 0.5;
			
			temp2->binary = 1;
			temp2->rel = 1;
		
			temp1->binary = 1;
			temp1->rel = 1;
            
			float second_call= factoring_irrel_C (0);
			
			temp1->binary = 0; //delete edge
			temp1->rel = 0.5;
            temp2->binary = 0;
			temp2->rel = 0.5;
			
			//cout << endl << "vertex: " << temp1->vertex;
			float first_call = factoring_irrel_C (1);
       
			
              temp1->binary = 1; //delete edge
		      temp1->rel = 0.5;
              temp2->binary = 1;
			  temp2->rel = 0.5;
			  if (flag && irr) // put back irrelevant edges
			  { for (int j=0; j < edge_counter ; j++)
			        { y = edge_irrel_1[j].v1;
			          r = edge_irrel_1[j].v2;
                      temp1 = &headnodes[y];//delete edge
		              while(temp1->vertex!= r)
			           {	
				         temp1 = temp1->next; // tranfer the address of 'temp->next' to 'temp'
			           }  
					  temp1->binary = 1;
					  temp1->rel = 0.5;
					  /*cout << endl << "======= inderting irrel ======";
					  /*cout << endl << "edge: (" << u <<", "<< r <<")";
                      /*cout << endl << "======= inderting irrel ======";*/
			  }
			  delete [ ] edge_irrel_2;
			  }
			return (0.5 * first_call + 0.5 * second_call);
		}
        
} 
	void create()
	//create function
	{   //cout << "create: " << endl;
		node *pre;
		node * nextn;
	    node *newnode;
		for (int i=0; i < n; i++)
		 for (int j=0; j < n; j ++)
		 {  
			infile >> A;
			//cout << "entry: " << A;
			if (i==j) A=A/2;
			weight [i][j] = A;
			if (A > 0 ) // exclude self-loops
			{ for (int k = 1; k <= A;  k++)
			   { newnode= new node;
                 newnode->vertex = j;
			     
                 if( headnodes[i].next == NULL )
			    {
			  	    newnode->next= NULL;
				
				    headnodes[i].next=newnode;
			   }
			else
			   {
				 pre= &headnodes[i];
				while( pre->next != NULL )
				{   
					 pre = pre->next;
				}
				newnode->next = NULL;
				pre->next = newnode;
				 }
		//ADJACENT NODES
		/*
			newnode= new node;
			newnode->vertex = i;
			if( headnodes[j].next == NULL )
			{
				newnode->next= NULL;
				headnodes[j].next=newnode;
			}
			else
			{
				pre= &headnodes[j];
				while( pre->next != NULL )
				{
					pre = pre->next;
				}
				newnode->next = NULL;
				pre->next = newnode;
			}
			*/
			} // for internal
			} // if
			} // for external
	
    
	}
	void DFS(int father, int v)
	// DFS function
	{  
	Visited [v]=true;
	bool adjtoa=true;
    node* adjnode=headnodes[v].next;
    while (adjnode) // visit all vertices adjacent to v
	{
		if (!Visited[adjnode->vertex])
		{//if adjacent vertex to v was not visited previously
                      DFS(v,adjnode->vertex);
		}
		else if (father !=adjnode->vertex) // if the vertex adjacent to v is not the father, we have a 
		{ 
			cycles=true;
        }
		adjnode = adjnode->next;
	}
	}
	void Bi_Connect(int father, int v)
	// DFS function
	{  
	  Visited [v]=true;
	  bool adjtoa=true;
	  edge b_e;
	  int b_min;
	  b_num [v] = b_num1; //tree edge 
	  b_low [v] = b_num [v]; b_num1 ++;
	
    node* adjnode=headnodes[v].next;
    while (adjnode) // visit all vertices adjacent to v
	{   if ((adjnode -> vertex !=father) && (b_num [adjnode -> vertex] < b_num [v])) // push edge to the stack
		//if ((b_num [adjnode -> vertex] < b_num [v])) // push edge to the stack
     	{  
		   b_e.v1 = v;
	       b_e.v2 = adjnode->vertex;
		         mystack.push (b_e);
		         mystack1.push (b_e); // auxiliary stack
		   
		   //cout << endl << "stack size: " << mystack.size(); //pushing it but not popping it.
	    }
	    if (!b_num[adjnode->vertex]) {
			Bi_Connect(v,adjnode->vertex);
		    if (b_low [v] > b_low [adjnode->vertex])
			    b_low [v] = b_low [adjnode->vertex];
			//cout << endl << "b-low of w= " << adjnode->vertex << "-" << b_low [adjnode->vertex] << "b_num v=" << v << "-" << b_num [v];
			      if (b_low [adjnode->vertex] >= b_num [v]) { 
					  
					//if (Exm_cmpnt (v, adjnode->vertex)) // PSEUDO-KNOT FOUND
					 // {
				          outfile << endl << "===================== New Block ================== \n" << endl;
				          do { 
						      
					           // delete an edge from the stack of the stack
				      	        b_e=mystack.top();
					            mystack.pop();
								for (int l=1; l <= weight [b_e.v1] [b_e.v2]; l++)
					                  outfile << "(" << b_e.v1 << "," << b_e.v2 <<") - ";
				                } while ( !((b_e.v1==adjnode->vertex && b_e.v2==v) || (b_e.v2==adjnode->vertex && b_e.v1==v))); 
					       outfile << endl;
						   if (Exm_cmpnt (v, adjnode->vertex)) {
	                            numpseudo_blocks ++;
                                outfile << endl << endl << " ---- this block represents a pseudoknot ---- " << endl ;
						   }
						   else {
						       numreg_blocks++;
						       outfile << endl << " ---- this block represents a regular-region ---- "<< endl;
						   }
					 // }
				  }
                 
		}
				  
		else if ( adjnode->vertex != father) 
		          {  if (b_num [adjnode->vertex] < b_low [v])
					 b_low [v]=b_num [adjnode->vertex];
		          }
		adjnode = adjnode->next;
	}
			
	}
	
	bool Exm_cmpnt (int v, int w)
    {    bool pseudoknot = false; 
	     edge b_e;
		
		         for (int i = 0 ; i < n ; i++) {
					 //T[i] = false;
					 d_cmpnt [i] = 0; 
				 }
		         do { 
						      
					           // delete an edge from the stack
				      	        b_e=mystack1.top();
					            mystack1.pop();
								for (int l=1; l <= weight [b_e.v1] [b_e.v2]; l++){ // parallel edges
					                  //outfile << "(" << b_e.v1 << "," << b_e.v2 <<") - ";
								     d_cmpnt [b_e.v1]++; 
								     d_cmpnt [b_e.v2]++;
		
								}
								if (d_cmpnt [b_e.v1] > 2 || d_cmpnt [b_e.v2] > 2) 
									{pseudoknot = true;//biconnected component has a vertex of degree 3ff
								    if (d_cmpnt [b_e.v1] > 2)
								     outfile << endl << "degree of " << b_e.v1 << " is " << d_cmpnt [b_e.v1] ;
									if (d_cmpnt [b_e.v2] > 2)
									 outfile << endl << "degree of " << b_e.v2<< " is " << d_cmpnt [b_e.v2] ;
								}
			    	 } while ( !((b_e.v1==w && b_e.v2==v) || (b_e.v2==w && b_e.v1==v))); 
					      
				return pseudoknot;
	}
	void setup_counter ()
	{ counter = 0;
	}
	long get_counter ()
	{ return counter;
	}
	void Dijkstra_B (int &val,edge &x1)//disjtra modified for factoring
		//Dijkstra function
       {   
		   bool all1s = true;
		   x1.v1= -1;
		   x1.v2= -1;
           int INFINITY = 1000;
		   int s =0,Minimum,u;
		   //int temp[40];			//temporary data structure which is needed to 
		   //T = new bool[n];
		   //L= new int[n];
		   //father = new int[n];
		   
		//INITIALIZING DATA STRUCTURES
		   for(int i =0;i<n;i++)
			{
			T[i] = false;			//T=V; all the vertices are eligible to be visited
			L[i]=INFINITY;			// at the beginning every vertex is at distance ??, from s
			temp[i] = INFINITY;		
			father[i] =-1;
			}
		   //WE ASSUME THE SOURCE IS 0 FOR EVERY GRAPH
		   L[s]=0;					// s is at distance 0 of itself
		   temp[s] =0;				// Temp has the same contents of L 
		  
		  // Let u be the vertex in T that has minimum L clearly at the beginning u=s 
		   for(int i = 0; i < n; i++)				//LOOP TROUGH ALL THE VERTICES OF THE GRAPH
		  {  
			  //cout<<endl<<"STEP "<<i<<":\n________ ";
			  Minimum = INFINITY;
		      for(int j = 0; j < n; j++)
			  {   if (!T[j]){
				  if( Minimum > temp[j])
					   {
						   Minimum = temp[j];			//finding minimum L[s]
						   u = j;
					   }
			   }
			  }
			  //temp[u] = INFINITY;				//Assigning INIFINITY to the data structure allready visited to find the next minimum L
			  //cout<<"\nU : "<<u;
			  w = &headnodes[u];
			  w=w->next;
			  //Assigning address of headnodes[u] to find adjacent nodes for u
			  //cout<<"\tW = "<<w->vertex;
			  while( w!=NULL )					//while w adjacent to u 
			  {
				  if(T[w->vertex]==false && w->binary)		// if w Exist in T, proceed
				  {   if (all1s && w->rel==0.5)// not all edges are 1
				       {all1s = false;
				          x1.v1=u;
						  x1.v2=w->vertex;
				       }
					  if (L[w->vertex]> L[u]+ weight[u][w->vertex])
					  { 
						  // if by reaching w from u has less weight
                                L[w->vertex]= L[u]+ weight[u][w->vertex]; // w is closer to s by using u;
								//cout<<"\nL[w]= L[u]+ weight(u,w)\n";
								//cout<<L[w->vertex]<<"  =   "<<L[u]<<"   +   "<<weight[u][w->vertex]<<endl;
								temp[w->vertex] = L[w->vertex];
                                father[w->vertex]=u;
								//cout<<"father[w] = "<<u<<endl;
						}
				  }
			
				w = w->next;   // tranfer the address of 'temp->next' to 'temp'
			  }
			  T[headnodes[u].vertex] = true;//Discard visited vertex u from T
			 // if (u == n-1)
              if (u == terminal)
				  break;	  
		   }
			  if (all1s)
				  //if (L[n-1] <= D)
				  if (L[terminal] <= D)
					  val = 1;
				  else val =0;
				  
			  else { // no all paths have value 1
				  //if ( L[n-1] > D )
				  if ( L[terminal] > D )
					   val = 0;
				  else val = 2;  //continue
			  }
			  
	}
	
	void display_val ()
	{   node *temp1;
		 for(int i = 0;i<n;i++)
		{
			temp1 = &headnodes[i];
			while( temp1->next!=NULL )
			{   if (temp1->next->binary)
			{
				cout << endl;
				cout<< "(" << i << "),(" << temp1->next->vertex<<")"; 
			}      
				temp1 = temp1->next; // tranfer the address of 'temp->next' to 'temp'
			}   
		}
	}
	int Dltirrlvnt(edge edge_irrel_1[], int &count_edges)
	{ 
	    node *temp1, *nextn;
		int u1;
	    
       /*for(int i = 0;i<n;i++)
		{
			temp1 = &headnodes[i];
			while( temp1->next!=NULL )
			{
				cout<< "(" << i << "),(" << temp1->next->vertex<<") binary: ";// show the data in the linked list
                cout<< temp1->next->binary<<"  " << endl;   
				
				temp1 = temp1->next; // tranfer the address of 'temp->next' to 'temp'
			}   
		}*/
        
		//int * Ls = new int [n];
		//int * Lt = new int [n];
		//edge *edge_irrel= new edge [e];
		int count = 0;
		Dijkstra_A(0, headnodes, Ls);
		//cout << "first call to Disjktras" << endl;
		Dijkstra_A(terminal, headnodes, Lt); 
		//cout << "second call to Disjktras" << endl;
			for (int i = 0; i < n; i++)
	       {
            nextn = &headnodes[i];
			nextn=nextn->next;
			
			while (nextn!=0)
			{  
			   u1=nextn->vertex;
			   if ((Ls[i]+Lt[u1] >=D) && (Ls[u1] + Lt[i]>=D)&&  nextn->binary==1 && nextn->rel==0.5){
				   nextn->binary=0;
                   nextn->rel = 0;
				   edge_irrel_1[count].v1=i;
                   edge_irrel_1[count].v2=u1;
				   count++;
				   
				   //cout << endl << "edge= (" << i << "," << u1 << "): " << Ls[i] + Lt[u1]; 
				   //cout << ", " << Ls[u1] + Lt[i];
				   //cout << "-----------------------" << endl;
			   }
			   nextn=nextn->next;
			}
		}
     
		count_edges = count;
		if (count > 0) return (1);
		else return (0);
				
	}	
int Dltirrlvnt_C(edge edge_irrel_1[], int &count_edges)
	{ 
	    node *temp1, *nextn;
		int u1;
	    
      
        
		
		int count = 0;
		
			for (int i = 0; i < n; i++)
	       {
            nextn = &headnodes[i];
			nextn=nextn->next;
			
			while (nextn!=0)
			{  
			   u1=nextn->vertex;
			   //Dijkstra_C(, headnodes, Ls);
			   if (nextn->binary==1 && nextn->rel==0.5 && i<u1){// possible irrel
				   //nextn->binary=0;
				   // find equivalent edge 
                   temp1 = &headnodes[u1];//delete edge
		           while(temp1->vertex!= i)
			   {	
				   temp1 = temp1->next; // tranfer the address of 'temp->next' to 'temp'
			   }  
		           // +++++++++ delete the edge (i,u) and (u,i) from G
			       temp1->binary = 0;
			       temp1->rel = 0.5;
                   nextn->binary=0;
				   nextn->rel=0.5;
				   Dijkstra_A(0,headnodes, Ls);
                   Dijkstra_A(terminal, headnodes, Lt);
				   if ((Ls[i]+Lt[u1] >=D) && (Ls[u1] + Lt[i]>=D))
				   {
                       edge_irrel_1[count].v1=i;
                       edge_irrel_1[count].v2=u1;
				       count++;
                       edge_irrel_1[count].v2=i;
                       edge_irrel_1[count].v1=u1;
				       count++;
				   } else// not irrelevant
				        {  temp1->binary = 1;
			               temp1->rel = 0.5;
                           nextn->binary=1;
				           nextn->rel=0.5; 
				  
				   }
			   }
			   nextn=nextn->next;
			} // while
			} // for
      //cout << endl << "********************updated list";
		//node *temp1;
		
		/*for(int i = 0;i<n;i++)
		{
			temp1 = &headnodes[i];
			while( temp1->next!=NULL )
			{
				//cout<< "(" << i << "),(" << temp1->next->vertex<<") binary: ";// show the data in the linked list
                //cout<< temp1->next->binary<<"  " << endl;
				temp1 = temp1->next; // tranfer the address of 'temp->next' to 'temp'
			}   
		}*/
		
		//cout<<"Number of vertices: "<<n<<endl;
        //cout<<"Number of edges: "<<e<<endl;
		count_edges = count;
		if (count > 0) return (1);
		else return (0);
				
}				
		void display()
	{
		node *temp1;
		outfile << endl << "==============================================";
		outfile << endl;
		
		for(int i = 0;i<n;i++)
		{
			temp1 = &headnodes[i];
			while( temp1!=NULL )
			{
				outfile << temp1->vertex<<" -> ";// show the data in the linked list
				temp1 = temp1->next;   // tranfer the address of 'temp->next' to 'temp'
			}
			outfile <<endl;
		}
		
		//cout<<"Number of vertices: "<<n<<endl;
        //cout<<"Number of edges: "<<e<<endl;
	
	}	
	
//;// end class
void Dijkstra_A (int a, node headnodes_A[],int L[])
{
		   int INFINITY = 1000;
		   int s =a,Minimum,u;
		   //int temp[30];			//temporary data structure which is needed to 
		   //T = new bool[n];
		   //L= new int[n];
		   //father = new int[n];
		   //int weight [n] [n];
		//INITIALIZING DATA STRUCTURES
		   for(int i =0;i<n;i++)
			{
			T[i] = false;			//T=V; all the vertices are eligible to be visited
			L[i]=INFINITY;			// at the beginning every vertex is at distance ??, from s
			temp[i] = INFINITY;		
			father[i] =-1;
			}
		   //WE ASSUME THE SOURCE IS 0 FOR EVERY GRAPH
		   L[s]=0;					// s is at distance 0 of itself
		   temp[s] =0;				// Temp has the same contents of L 
		  
		  // Let u be the vertex in T that has minimum L clearly at the beginning u=s 
		   for(int i = 0; i < n; i++)				//LOOP TROUGH ALL THE VERTICES OF THE GRAPH
		  {  
			  //cout<<endl<<"STEP "<<i<<":\n________ ";
			  Minimum = INFINITY;
		      for(int j = 0; j < n; j++)
			  {   if (T[j]==0){
				  if( Minimum > temp[j])
					   {
						   Minimum = temp[j];			//finding minimum L[s]
						   u = j;
					   }
			      }
			  }
			  //temp[u] = INFINITY;				//Assigning INIFINITY to the data structure allready visited to find the next minimum L
			  //cout<<"\nU : "<<u;
			  w = &headnodes_A[u];
			  w = w->next;
			  //Assigning address of headnodes[u] to find adjacent nodes for u
			  //cout<<"\tW = "<<w->vertex;
			  while( w!=NULL )					//while w adjacent to u 
			  {
				  if(T[w->vertex]==false && w->binary)		// if w Exist in T, proceed
				  {
					  if (L[w->vertex]> L[u]+ weight[u][w->vertex])
					  { 
						  // if by reaching w from u has less weight
                                L[w->vertex]= L[u]+ weight[u][w->vertex]; // w is closer to s by using u;
								//cout<<"\nL[w]= L[u]+ weight(u,w)\n";
								//cout<<L[w->vertex]<<"  =   "<<L[u]<<"   +   "<<weight[u][w->vertex]<<endl;
								temp[w->vertex] = L[w->vertex];
                                father[w->vertex]=u;
								//cout<<"father[w] = "<<u<<endl;
						}
					}
			
				w = w->next;   // tranfer the address of 'temp->next' to 'temp'
			  }
			  T[u] = true;				//Discard visited vertex u from T
		
		   }
}
};// end class
int main(int argc, char ** argv)
{
    string input, output;
    //cout<<argv[0]<<endl;
    if(strcmp(argv[1],"-input")==0)
    {
        input=argv[2];
    }
    
    
    if(strcmp(argv[3],"-len")==0)
    {
        nods=atoi(argv[4]);
    }
    
    if(strcmp(argv[5],"-output")==0)
    {
        output=argv[6];
    }
    infile.open(input);
	outfile.open(output);
    
    //infile.open("/Users/cs4367/Documents/LouisLab/4OJI_adj.txt");
	//outfile.open("/Users/cs4367/Documents/LouisLab/4OJI_adj_out.txt");
    
    if (infile.is_open()){
        
	clock_t start;
	clock_t end;
	float facto_time;
	float facto_irr_time;
	float facto_irr_e_time;
	string motif_number;
	//cout<< endl << "Please enter the number of vertices: ";
	//cin >>nods;
    
	
	//cout<< endl << "Please enter motif number: ";
	//cin >> motif_number;
    
	edgs=0;
	dim=0;
	term = 0;
    Visited=new bool [nods];
	int b_counter = 0;
	
    int Nmbrcmpnts=0; //we initialize the counter for the number of components
    
    if  (outfile.is_open()){
	
	//outfile << endl << "--------------------- Motif :" << motif_number << " -----------------------------" << endl;
    }
        else{
            cout<<"Cannot open outputfile!"<<endl;
            exit(0);
        }
    //outfile << " +++++++++++++++ number of vertices: " << nods << "    ====================" << endl;
	Graph G(nods,edgs,dim,term);
	while ( !infile.eof () && b_counter < 1)
	{G.init();
	    //cout << endl << "ok";
	    b_counter ++;
		
	    G.create ();
		 //cout << endl << "ok2";
		//G.display();
		//cout << endl << "ok3";
		
        G.Bi_Connect(-1, 0);
		//outfile << endl << "----------- Summary information for Motif :" << motif_number << " --------------------------------" << endl;
		outfile << "----------- Total number of blocks: " << numreg_blocks + numpseudo_blocks << endl;
		outfile << "----------- number of PK blocks: " << numpseudo_blocks << endl;
		outfile << "----------- number of regular blocks : " << numreg_blocks << endl;
		outfile << "-------------------------------------------------------------------------------------" << endl;
	} 
	
    infile.close();
	outfile.close();
    //system("pause");
    }
    else{
        cout<<"Cannot open inputfile!"<<endl;
        exit(0);
    }
  return 0;
 
}
