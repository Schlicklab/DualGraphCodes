//Petingi 08/06/14
//to detect pseudo-graph in dual graphs
//collaboraton project with Tamar Schlick and Namhee Kim (NYU)
//we are reading AdJ Matrix thus we convert to Linked List DS
//Correction Swati found a mistake in the code (the last parallel edges were counted twice). 11/11/17

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <ctime>
#include <stack>
#include <string.h>
#include <stdlib.h>
#include <map>
#include <set>

using namespace std;


//node
struct edge
{
    int v1;
    int v2;
};
struct node
{
    int vertex;
//	int binary;
//	float rel;
    node * next;
	
};

struct subgraph_node // S.J. 03/29/2018 - to keep track of subgraphs
{
    int num_blocks; // 04/06/2018 - to keep track of number of basic blocks this graph is made up of
    vector<bool>vertices;
    bool pknot;
    vector<int>artPointsCount; // 04/03/2018 - to store the articulation points in this subgraph and how many times they occur - will be updated as we join graphs
};

int nods = 0;
int motif_counter = 0;
int numreg_blocks = 0;
int numpseudo_blocks = 0;
int edgs = 0;
int dim = 0;
int iterat = 0;
int term = 0;
stack<edge> mystack;
stack<edge> mystack1;
bool * Visited;
int *d_cmpnt; 
node * glovalnode;
int weight[1000][1000];
ifstream infile;
ofstream outfile;
bool cycles=false;// determine if a graph has cycles.
vector<subgraph_node>Blocks; // S.J. 03/29/2018
vector<int>artPoints; // 04/03/2018 - to store articulation points and how many times they occur
map<int,vector<subgraph_node> > artPointsBlocks; // to store basic blocks that correspond to each subgraph
bool all_subgraphs=false; // S.J. 03/30/2018 - flag to see if we want to calculate all subgraphs or not
int maxVertex=15; // S.J. 10/18/2018 - maximum vertices for which we want to calculate subgraphs


// S.J. 04/03/2018 - function to compare two subgraphs and see if they are same or not
//returns true when two subgraphs are not the same, so that it is inserted in the set Subgraphs
//returns false when two subgraphs are the same, so that it is not inserted
// 04/05/2018 - changed the code to provide a strict ordering of set elements - the graphs with more vertices are "greater"; if same number of vertices. then graphs with the vertex with the smallest number is "smaller"
bool comp_subgraphs(subgraph_node node1, subgraph_node node2){
    
    /*cout << "Subgraph:";
    for(int i=0; i<nods; i++){
        
        if(node1.vertices[i])
            cout << "\t" << i;
    }
    cout << "end" << endl;
    
    cout << "Subgraph:";
    for(int i=0; i<nods; i++){
        
        if(node2.vertices[i])
            cout << "\t" << i;
    }
    cout << "end" << endl;
    
    for(int i=0;i<nods;i++){
        
        // the subgraphs are not the same as one if true and other is false for this vertex
        if(node1.vertices[i] != node2.vertices[i]){
            cout << "I return true" << endl;
            return true;
        }
    }
    
    cout << "I return false" << endl;
    return false; // the two subgraphs are the same*/
    
    //04/05/2018 - new code for comparison as a strict ordering is needed.
    int num1=0,num2=0,sum1=0,sum2=0,diff=-1; // 04/14/2013 - to change the condition for strict lower if same number of vertices
    
    for(int i=0;i<nods;i++){
        
        if(node1.vertices[i]){
            num1++;
            //sum1+=i;
        }
        if(node2.vertices[i]){
            num2++;
            //sum2+=i;
        }
        
        if(node1.vertices[i] != node2.vertices[i] && diff == -1){ // 04/19/2018 - change to make sure we only get the first vertex they are different
            diff=i; // the first vertex where the two subgraphs are different
        }
    }
    if(num1!=num2)
        return num1<num2;
    else {
        //return sum1<sum2;
        if(diff !=-1) // 04/19/2018 - when two graphs are not the same - the graph with the smaller vertex number is smaller (prence of a vertex is 1 and absence is 0, therefore > sign instead of <)
            return node1.vertices[diff]>node2.vertices[diff];
        else // two graphs are the same
            return node1.vertices[0]<node2.vertices[0];
    }
    
}
bool(*fn_pt)(subgraph_node,subgraph_node) = comp_subgraphs;
set<subgraph_node,bool(*)(subgraph_node,subgraph_node)> Subgraphs (fn_pt);

//S.J. 03/30/2018 - insert a new block into the list
/*void insert_subgraph(subgraph_node newnode){
    
    vector<subgraph_node>::iterator it;
    
    for(it=Blocks.begin();it!=Blocks.end();it++){
        
        for(int i=0;i<nods;i++){
            if(newnode.vertices[i] && it->vertices[i]){ // there is a common vertex, insert this in between
                
                Blocks.insert(it+1,newnode);
                return;
            }
        }
    }
    
    // if the code comes here, then there are no common vertices with current blocks/or this is the first block, add the newnode in the end
    Blocks.push_back(newnode);
    return;
}*/

//S.J. 04/03/2018 - function to setup the articulation points and corresponding graphs
void setup_artPoints(){
    
    int num_vertices=0; // 04/19/2018 to count number of vertices in this block
    
    //find articulation points in each subgraph, and update the artPointsCount in each, and all each art point with its corresponding blocks
    artPoints[0]--; // because 0 is always returned by the function
    for(int i=0;i<Blocks.size();i++){
        
        num_vertices=0;
        Blocks[i].artPointsCount.resize(nods); // initialize it with zero
        for(int j=0;j<nods;j++){
            
            if(artPoints[j]!=0 && Blocks[i].vertices[j]) // if any of the articulation points are in this subgraph, update the Subgraph
                Blocks[i].artPointsCount[j] = artPoints[j];
        }
        
        // has to be done after the above loop as the Ssubgraphs need to be updated fully
        for(int j=0;j<nods;j++){
            
            if(Blocks[i].vertices[j]){ // if this vertex exists in this subgraph - 04/19/2018
                
                num_vertices++;
                if(artPoints[j]!=0) // if this vertx is an articulation point, push this graph to its articulation point map
                    artPointsBlocks[j].push_back(Blocks[i]);
            }
        }
        
        if(num_vertices >=2 && num_vertices <=maxVertex) // 04/19/2018 - only add subgraphs if they are between 2-9 vertices, 10/18/2018 - changed it to the variable maxVertex
            Subgraphs.insert(Blocks[i]);
    }
    
    // print the basic blocks allocated to all articulation points
    /*for(int i=0;i<nods;i++){
     
        if(artPoints[i] !=0){ // this is an articulation point
     
            cout << "Articulation Point: " << i << endl;
            for(int j=0;j<artPointsBlocks[i].size();j++){
     
                cout << "Subgraph:";
                for(int k=0; k<nods; k++){
     
                    if(artPointsBlocks[i][j].vertices[k])
                        cout << "\t" << k;
                }
                cout << "\tPseudoknot:\t";
                if(artPointsBlocks[i][j].pknot)
                    cout << "Yes" << endl;
                else
                cout << "No" << endl;
            }
        }
     }*/
    
    // insert the basic blocks into the Subgraphs set
    //for(int i=0;i<Blocks.size();i++){
        
        //Subgraphs.insert(Blocks[i]);
        //cout << "Inserting block " << i << endl;
    //}

}

// S.J. 03/30/2018 - calculate and print all possible subgraphs and add it to the end of the list
//void calc_possible_subgraphs(){
void calc_possible_subgraphs(subgraph_node graph){ // 04/03/2018 - for the recursive function
    
    //Earlier code - valid for only linear graphs (graphs joined tother linearly, one after another)
    /*int total_blocks=0;
    subgraph_node newnode;
    vector<subgraph_node> higher_subgraphs;
    vector<subgraph_node>::iterator start, temp, temp2;
    
    total_blocks = numreg_blocks + numpseudo_blocks; // total number of blocks

    for(int i=2;i<=total_blocks;i++){ // for loop to increase the number of blocks in the new subgraph
        
        start=Blocks.begin();
        for(int j=0;j<total_blocks-i+1;j++){ // starting point for the calculation of the subgraphs
            
            //clearing the previous contents of the newnode
            newnode.vertices.clear();
            newnode.vertices.resize(nods);
            newnode.pknot=false;
            
            temp=start;
            for(int k=1;k<=i;k++){ // for loop to include number of blocks
             
                for(int m=0;m<nods;m++){ // add the vertices of the blocks that are part of the this subgraph
                    
                    if(temp->vertices[m])
                        newnode.vertices[m]=true;
                }
                if(temp->pknot) // record if any of the blocks of this subgraph contain pseudoknots
                    newnode.pknot=true;
                
                temp++;
            }
            //insert new graph at the begining of the list
            higher_subgraphs.insert(higher_subgraphs.begin(),newnode);
            start++;
        }
    }
    if(higher_subgraphs.size() != 0) // if any new subgraphs are created
        Blocks.insert(Blocks.begin(),higher_subgraphs.begin(),higher_subgraphs.end());  // making one combined list
    */
    
    //04/03/2018 - new code for all possible subgraphs
    subgraph_node newsubgraph;
    int num_vertices=0; // 04/14/2018 - to count the number of vertices in the subgraph
    
    for(int i=0;i<nods;i++){
     
        if(graph.artPointsCount[i] != 0){ // for every articulation point of this subgraph
            
            //cout << "Articulation Point: " << i << endl;
            for(int j=0; j<artPointsBlocks[i].size(); j++){ // for every block associated with this articulation point
                
                //cout << "Art points block number: " << j << endl;
                //resetting the subgraph structure
                newsubgraph.vertices.clear();
                newsubgraph.vertices.resize(nods);
                newsubgraph.artPointsCount.clear();
                newsubgraph.artPointsCount.resize(nods);
                newsubgraph.pknot=false;
                num_vertices=0;
                
                for(int k=0;k<nods;k++){ //combine the vertices and artPointsCount of the two graphs
                    
                    if(graph.vertices[k] || artPointsBlocks[i][j].vertices[k]){ // if either graph has this vertex
                        newsubgraph.vertices[k]=true;
                        num_vertices++; //04/14/2018
                    }
                    
                    if(k!=i) // if this vertex is not the current articulation point
                        newsubgraph.artPointsCount[k]=graph.artPointsCount[k]+artPointsBlocks[i][j].artPointsCount[k];
                    else // the count of this articulation point should decrease by one
                        newsubgraph.artPointsCount[k]=graph.artPointsCount[k]-1;
                }
                newsubgraph.num_blocks=graph.num_blocks+1; // 04/06/2018 - as the new graphs are always created by combining the current graph with a block
                if(graph.pknot || artPointsBlocks[i][j].pknot)
                    newsubgraph.pknot=true;
                
                if(num_vertices > maxVertex) // 04/14/2018 if the number of vertices in greater than 9, then just ignore it, 10/18/2018 - changed it to variable maxVertex
                    continue;
                
                //cout << "NewSubgraph:";
                //for(int k=0; k<nods; k++){
                //    if(newsubgraph.vertices[k])
                //        cout << "\t" << k;
                //}
                //cout << "end" << endl;
                
                //if this subgraph does not exist in the set already
                if(Subgraphs.find(newsubgraph) != Subgraphs.end()){
                    //cout << "Subgraph already found" << endl;
                    continue;// nothing to do if it already exists
                }
                else{
                
                    Subgraphs.insert(newsubgraph); // insert this in the Subgraphs
                    //cout << "Subgraph inserted: Calling recursion" << endl;
                    if(num_vertices<maxVertex) // only call the recursion if the size of the subgraph is < 9, 10/18/2018 - changed it to maxVertex
                        calc_possible_subgraphs(newsubgraph); // recurse on the new subgraph
                }
            }
        }
    }
    
    return;
}

//S.J. 03/30/2018 - to print the subgraphs
void print_subgraphs(string outfile_all){
    
    ofstream output;
    output.open(outfile_all);
    //vector<subgraph_node>::iterator temp;
    set<subgraph_node,bool(*)(subgraph_node,subgraph_node)>:: iterator temp;
    int numreg=0, numpk=0;

    //for(temp=Blocks.begin();temp!=Blocks.end();temp++){
    for(temp=Subgraphs.begin();temp!=Subgraphs.end();temp++){
        
        /*output << "Subgraph:";
        for(int i=0; i<nods; i++){
            
            if(temp->vertices[i])
                output << "\t" << i;
        }
        output << "\tPseudoknot:\t";
        if(temp->pknot)
            output << "Yes\t";
        else
            output << "No\t";
        
        output << "Num Blocks:\t" << temp->num_blocks << endl;*/
        
        // 04/14/2018 - printing the subgraphs like blocks
        output << endl << "===================== New Subgraph ================== \n" << endl;
        for(int i=0;i<nods;i++){
            
            if(temp->vertices[i]){ // if this vertex exists in this graph
                
                for(int j=i+1;j<nods;j++){
                    
                    if(temp->vertices[j]){ // if this vertex exists in this graph
                
                        for(int k=1;k<=weight[i][j];k++) // print the number of edges
                            output << "(" << i << "," << j << ") - ";
                    }
                }
            }
        }
        if(temp->pknot){
            output << endl << endl << " ---- this subgraph contains a pseudoknot ---- " << endl ;
            numpk++;
        }
        else{
            output << endl << endl << " ---- this subgraph represents a regular-region ---- "<< endl;
            numreg++;
        }
    }
    
    //04/19/2018 - to make it consistent with
    output << "----------- Total number of subgraphs between 2-15 vertices: " << numreg + numpk << endl;
    output << "----------- number of PK containing subgraphs: " << numpk << endl;
    output << "----------- number of regular subgraphs : " << numreg << endl;
    output << "-------------------------------------------------------------------------------------" << endl;

    
    output.close();
}

//S.J. 03/30/2018 - to parse the command line arguments
void parse_command_line(int argc, char ** argv, string & input, string & output, string & outfile_all){
    
    for(int i=1;i<argc;i++){
    
        if(strcmp(argv[i],"-input")==0)
        {
            input=argv[i+1];
        }
        else if(strcmp(argv[i],"-len")==0)
        {
            nods=atoi(argv[i+1]);
        }
        else if(strcmp(argv[i],"-output")==0)
        {
            output=argv[i+1];
        }
        else if(strcmp(argv[i],"-all")==0)
        {
            all_subgraphs=true;
            outfile_all=argv[i+1];
        }
    }
}

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
    
	Graph (int nods, int edgs, int dim, int term) // construtor
    {
		//infile.open("C:/users/petingi/desktop/RnaDv/V7adj.txt");
        terminal=term;
        //counter;
		n=nods;
		b_num1 = 0;
		b_low1 = 0;
		b_low = new int [n];
		b_num = new int [n];
		
		e = edgs;
		D = dim;
		counter=0;
        headnodes= new node [n];
        //int temp[200];			//temporary data structure which is needed to
		T = new bool[n];
		L = new int[n];
		father = new int[n];
		d_cmpnt = new int [n];
        temp = new int[n];
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
    
    void init()// construtor
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
	
    ~Graph()
	{
        delete [] headnodes;
		//int temp[200];			//temporary data structure which is needed to 
		delete [] T;
		delete [] L;
		delete [] father;
        delete [] edge_irrel;
        delete [] Ls;
        delete [] Lt;
    }
    
	void Del( )
	{
        delete [] headnodes;
	    //delete [] headnodes_A;
		delete [] T;
		delete [] L;
		delete [] father;
        delete [] edge_irrel;
        delete [] Ls;
        delete [] Lt;
        //delete [][] weight;
    }
	
/*    void setdiam (int diameter)
    {
        D=diameter;
    }
    
	float factoring (int flag)
	{   
		int value;
		++counter;
		
	    Dijkstra_B (value, x);
		//cout << endl << "value: " << value << " " << x.v1 << " " << x.v2;
		if (value == 1) return 1;
		else if (value==0) return 0;
		else{
            // value = 2
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
		else {
            
            if (flag)  {
                
                edge_irrel_2=new edge [2*e];
			    irr=Dltirrlvnt(edge_irrel_2, edge_counter);
				if (!irr) delete [ ] edge_irrel_2;
				if (irr) {
                    
                    edge_irrel_1 = edge_irrel_2;
					Dijkstra_B (value, x1);
                    //cout << endl << "+++++++++++++++++++++++++++";
                    //cout << endl << "cuando hay irrel edge " << x1.v1 << "-" << x1.v2;
                    if (value < 2){
                       
                        for (int j=0; j < edge_counter ; j++){
                            
                            y = edge_irrel_1[j].v1;
                            r = edge_irrel_1[j].v2;
                            temp1 = &headnodes[y];//delete edge

		                    while(temp1->vertex!= r){
                                
                                temp1 = temp1->next; // tranfer the address of 'temp->next' to 'temp'
			                }  
                            temp1->binary = 1;
					        temp1->rel = 0.5;
					        //cout << endl << "======= inderting irrel ======";
					        //cout << endl << "edge irrelevante: (" << y <<", "<< r <<")";
                            //cout << endl << "======= inderting irrel ======";
                        }// for
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

            while(temp1->vertex!= x1.v2){
				temp1 = temp1->next; // tranfer the address of 'temp->next' to 'temp'
            }
		
            temp1->binary = 1;
            temp1->rel = 1;

            temp2 = &headnodes[x1.v2];//delete edge

            while(temp2->vertex!= x1.v1){
				temp2 = temp2->next; // tranfer the address of 'temp->next' to 'temp'
            }  //delete edge
			//temp2->rel = 0.5;
			
		
			temp2->binary = 1;
			temp2->rel = 1;
			//cout << endl << "primero: " << temp2->vertex;
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
				  
				  for (int j=0; j < edge_counter ; j++){
                      
                      y = edge_irrel_1[j].v1;
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
		else {
            
            if (flag) {
                
			    edge_irrel_2 = new edge [2*e];
	            irr=Dltirrlvnt_C(edge_irrel_2,edge_counter);
				if (!irr) delete [] edge_irrel_2;
                if (irr){
                    
                    edge_irrel_1= edge_irrel_2;
                    Dijkstra_B (value, x1);
                    if (value < 2) {
                        //cout << endl << "retorno counter: " << edge_irrel;
                        for (int j=0; j < edge_counter ; j++){
                            
                            y = edge_irrel_1[j].v1;
			                r = edge_irrel_1[j].v2;
                            temp1 = &headnodes[y];//delete edge

		                    while(temp1->vertex!= r){
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
            {
                for (int j=0; j < edge_counter ; j++){
                    
                    y = edge_irrel_1[j].v1;
			        r = edge_irrel_1[j].v2;
                    temp1 = &headnodes[y];//delete edge

		            while(temp1->vertex!= r){
                        
				         temp1 = temp1->next; // tranfer the address of 'temp->next' to 'temp'
                    }
					temp1->binary = 1;
					temp1->rel = 0.5;
					//cout << endl << "======= inderting irrel ======";
					//cout << endl << "edge: (" << u <<", "<< r <<")";
                    //cout << endl << "======= inderting irrel ======";
                }
                delete [ ] edge_irrel_2;
            }
			return (0.5 * first_call + 0.5 * second_call);
		}
    }
 */
	void create() //create function
	{
        //cout << "create: " << endl;
		node *pre;
		node * nextn;
	    node *newnode;
        
        for (int i=0; i < n; i++){
            
            for (int j=0; j < n; j++)
            {
                infile >> A;
                //cout << "entry: " << A;

                if (i==j) A=A/2;
                weight [i][j] = A;
                //cout << endl << "aqui";
                if (A > 0 ) // exclude self-loops // S.J. 03/28/2018 - it does not exclude self loops
                {
                    for (int k = 1; k <= 1;  k++){
                        
                        newnode= new node;
                        newnode->vertex = j;
                        //cout << endl << "here I am";
                        if( headnodes[i].next == NULL ){
                            
                            newnode->next= NULL;
                            headnodes[i].next=newnode;
                        }
                        else{
				 
                            pre= &headnodes[i];
                            while( pre->next != NULL ){
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
        } // first for loop
        
        //S.J. 03/28/2018 - printing the linked list
        /*node * temp;
        cout << "Pringing the linked list" << endl;
        for(int i=0; i<n; i++){
            
            temp = &headnodes[i];
            while(temp != NULL){
                
                cout << temp->vertex << "\t";
                temp=temp->next;
            }
            cout << endl;
        }
        
        cout << "Printing the weights" << endl;
        for(int i=0; i<n; i++){
            
            for(int j=0; j<n; j++){
                
                cout << weight[i][j] << "\t";
            }
            cout << endl;
        }*/
	}
    
/*	void DFS(int father, int v)// DFS function
	{  
        Visited [v]=true;
        bool adjtoa=true;

        node* adjnode=headnodes[v].next;
        while (adjnode) // visit all vertices adjacent to v
        {
            if (!Visited[adjnode->vertex]){//if adjacent vertex to v was not visited previously
                
                DFS(v,adjnode->vertex);
            }
            else if (father !=adjnode->vertex) // if the vertex adjacent to v is not the father, we have a
            {
                cycles=true;
            }
            adjnode = adjnode->next;

        }
	}
*/
	void Bi_Connect(int father, int v) // DFS function
	{  
        Visited [v]=true;
        bool adjtoa=true;
        edge b_e;
        int b_min;
        b_num [v] = b_num1; //tree edge
        b_low [v] = b_num [v]; b_num1 ++;
	
        subgraph_node new_subgraph; // S.J. 03/30/2018
        
        node* adjnode=headnodes[v].next;
        while (adjnode) // visit all vertices adjacent to v
        {
            if ((adjnode -> vertex !=father) && (b_num [adjnode -> vertex] < b_num [v])) // push edge to the stack
            //if ((b_num [adjnode -> vertex] < b_num [v])) // push edge to the stack
            {
                b_e.v1 = v;
                b_e.v2 = adjnode->vertex;
                //cout << endl << "edge: (" << b_e.v1 << "," << b_e.v2 << ")";
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
                    //cout << "Vertex: " << v << endl;
                    outfile << endl << "===================== New Block ================== \n" << endl;
                    
                    // S.J. 03/29/2018 - create a new subgraph node
                    if(all_subgraphs){ // if all subgraphs have to be calculated, then only store the blocks and the articulation points
                        new_subgraph.vertices.clear();
                        new_subgraph.vertices.resize(nods);
                        new_subgraph.num_blocks=1; // 04/06/2018
                        artPoints[v]++; // 04/03/2018 - to store the number of times a articulation point comes, that determines how many subgraphs to merge at this point
                        
                    }
                    
                    do {
					     // delete an edge from the stack of the stack
                        b_e=mystack.top();
					    mystack.pop();
						//cout << endl << b_e.v1 << "'" << b_e.v2 << "  - weight: " << weight[b_e.v1][b_e.v2];
                        
                        // S.J. 03/29/2018 - adding the vertices that correspond to this subgraph to the node
                        if(all_subgraphs){
                            new_subgraph.vertices[b_e.v1]=true;
                            new_subgraph.vertices[b_e.v2]=true;
                        }
                        
                        for (int l=1; l <= weight [b_e.v1] [b_e.v2]; l++)
                            outfile << "(" << b_e.v1 << "," << b_e.v2 <<") - ";
                    } while ( !((b_e.v1==adjnode->vertex && b_e.v2==v) || (b_e.v2==adjnode->vertex && b_e.v1==v)));
                    
                    outfile << endl;
					if (Exm_cmpnt (v, adjnode->vertex)) {
                        numpseudo_blocks ++;
                        outfile << endl << endl << " ---- this block represents a pseudoknot ---- " << endl ;
                        
                        if(all_subgraphs)
                            new_subgraph.pknot=true;// S.J. 03/29/2018

					}
					else {
                        numreg_blocks++;
						outfile << endl << " ---- this block represents a regular-region ---- "<< endl;
                        
                        if(all_subgraphs)
                            new_subgraph.pknot=false; // S.J. 03/29/2018
					}
					 // }
                    if(all_subgraphs) // S.J. 03/20/2018 - insert the new subgraph into the list
                        //insert_subgraph(new_subgraph);
                        Blocks.push_back(new_subgraph); // 04/03/2018 - as for the new algo we don't need it inserted in a specific place
                }
            }
			else if ( adjnode->vertex != father){
                
                if (b_num [adjnode->vertex] < b_low [v])
					 b_low [v]=b_num [adjnode->vertex];
            }
            adjnode = adjnode->next;
        }
	}
	
	bool Exm_cmpnt (int v, int w)
    {
        bool pseudoknot = false;
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
			if (d_cmpnt [b_e.v1] > 2 || d_cmpnt [b_e.v2] > 2){
                
                pseudoknot = true;//biconnected component has a vertex of degree 3ff
				if (d_cmpnt [b_e.v1] > 2)
                    outfile << endl << "degree of " << b_e.v1 << " is " << d_cmpnt [b_e.v1] ;
				if (d_cmpnt [b_e.v2] > 2)
                    outfile << endl << "degree of " << b_e.v2<< " is " << d_cmpnt [b_e.v2] ;
            }
		} while ( !((b_e.v1==w && b_e.v2==v) || (b_e.v2==w && b_e.v1==v)));
					      
		return pseudoknot;
	}

	void setup_counter ()
	{
        counter = 0;
	}
    
	long get_counter ()
	{
        return counter;
	}
    
/*	void Dijkstra_B (int &val,edge &x1)//disjtra modified for factoring //Dijkstra function
    {
        bool all1s = true;
		x1.v1= -1;
		x1.v2= -1;
        int INF = 1000;
		int s =0,Minimum,u;
		//int temp[40];			//temporary data structure which is needed to
		//T = new bool[n];
		//L= new int[n];
		//father = new int[n];
		   
		//INITIALIZING DATA STRUCTURES
		for(int i =0;i<n;i++){
            
			T[i] = false;			//T=V; all the vertices are eligible to be visited
			L[i]=INF;			// at the beginning every vertex is at distance , from s
			temp[i] = INF;		
			father[i] =-1;
        }
		//WE ASSUME THE SOURCE IS 0 FOR EVERY GRAPH
		L[s]=0;					// s is at distance 0 of itself
		temp[s] =0;				// Temp has the same contents of L
		  
        // Let u be the vertex in T that has minimum L clearly at the beginning u=s
		for(int i = 0; i < n; i++)				//LOOP TROUGH ALL THE VERTICES OF THE GRAPH
        {
            //cout<<endl<<"STEP "<<i<<":\n________ ";
			Minimum = INF;
		    for(int j = 0; j < n; j++){
                if (!T[j]){
                    if( Minimum > temp[j]){
                        
                        Minimum = temp[j];			//finding minimum L[s]
						u = j;
                    }
                }
            }
			//temp[u] = INFINITY;	//Assigning INIFINITY to the data structure allready visited to find the next minimum L
            //cout<<"\nU : "<<u;
			w = &headnodes[u];
			w=w->next;
			//Assigning address of headnodes[u] to find adjacent nodes for u
			//cout<<"\tW = "<<w->vertex;
			while( w!=NULL )					//while w adjacent to u
            {
                if(T[w->vertex]==false && w->binary)		// if w Exist in T, proceed
				{
                    if (all1s && w->rel==0.5)// not all edges are 1
                    {
                        all1s = false;
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

        if (all1s){
		//if (L[n-1] <= D)
            if (L[terminal] <= D)
                val = 1;
            else
                val =0;
        }
        else { // no all paths have value 1
            //if ( L[n-1] > D )
			if ( L[terminal] > D )
                val = 0;
            else
                val = 2;  //continue
        }
	}
	
	void display_val ()
	{
        node *temp1;
		for(int i = 0;i<n;i++)
		{
			temp1 = &headnodes[i];
			while( temp1->next!=NULL ){
            
                if (temp1->next->binary){
				
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
	    
        //for(int i = 0;i<n;i++)
		//{
		//	temp1 = &headnodes[i];
		//	while( temp1->next!=NULL )
		//	{
		//		cout<< "(" << i << "),(" << temp1->next->vertex<<") binary: ";// show the data in the linked list
        //      cout<< temp1->next->binary<<"  " << endl;
		//
		//		temp1 = temp1->next; // tranfer the address of 'temp->next' to 'temp'
		//	}
		//}
        
		//int * Ls = new int [n];
		//int * Lt = new int [n];
		//edge *edge_irrel= new edge [e];
		int count = 0;
		Dijkstra_A(0, headnodes, Ls);
		//cout << "first call to Disjktras" << endl;

		Dijkstra_A(terminal, headnodes, Lt); 
		//cout << "second call to Disjktras" << endl;

        for (int i = 0; i < n; i++){
            
            nextn = &headnodes[i];
			nextn=nextn->next;
			
			while (nextn!=0){
                
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
		if (count > 0)
            return (1);
		else
            return (0);
    }

    int Dltirrlvnt_C(edge edge_irrel_1[], int &count_edges)
	{ 
	    node *temp1, *nextn;
		int u1;
	    int count = 0;
		
        for (int i = 0; i < n; i++){
            
            nextn = &headnodes[i];
			nextn=nextn->next;
			
			while (nextn!=0){
                
			   u1=nextn->vertex;
			   //Dijkstra_C(, headnodes, Ls);
			   if (nextn->binary==1 && nextn->rel==0.5 && i<u1){// possible irrel
				   //nextn->binary=0;
				   // find equivalent edge 
                   temp1 = &headnodes[u1];//delete edge

		           while(temp1->vertex!= i){
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

				   }
                   else// not irrelevant
				   {
                       temp1->binary = 1;
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
		
		//for(int i = 0;i<n;i++)
		//{
		//	temp1 = &headnodes[i];
		//	while( temp1->next!=NULL )
		//	{
				//cout<< "(" << i << "),(" << temp1->next->vertex<<") binary: ";// show the data in the linked list
                //cout<< temp1->next->binary<<"  " << endl;
		//		temp1 = temp1->next; // tranfer the address of 'temp->next' to 'temp'
		//	}
		//}
		
		//cout<<"Number of vertices: "<<n<<endl;
        //cout<<"Number of edges: "<<e<<endl;
		count_edges = count;
		if (count > 0)
            return (1);
		else
            return (0);
    }
		
    void display()
	{
		node *temp1;
		cout << endl << "==============================================";
		cout << endl;
		
		for(int i = 0;i<n;i++)
		{
			temp1 = &headnodes[i];
			while( temp1!=NULL )
			{
				cout << temp1->vertex<<" -> ";// show the data in the linked list
				temp1 = temp1->next;   // tranfer the address of 'temp->next' to 'temp'
			}
			cout <<endl;
		}
		
		//cout<<"Number of vertices: "<<n<<endl;
        //cout<<"Number of edges: "<<e<<endl;
	}

	void Dijkstra_A (int a, node headnodes_A[],int L[])
    {
        int INF= 1000;
		int s =a,Minimum,u;
		//int temp[30];			//temporary data structure which is needed to
		//T = new bool[n];
		//L= new int[n];
		//father = new int[n];
		//int weight [n] [n];
		//INITIALIZING DATA STRUCTURES
		for(int i =0;i<n;i++){
            
			T[i] = false;			//T=V; all the vertices are eligible to be visited
			L[i]=INF;			// at the beginning every vertex is at distance , from s
			temp[i] = INF;		
			father[i] =-1;
        }
		//WE ASSUME THE SOURCE IS 0 FOR EVERY GRAPH
		L[s]=0;					// s is at distance 0 of itself
		temp[s] =0;				// Temp has the same contents of L
		  
		// Let u be the vertex in T that has minimum L clearly at the beginning u=s
		for(int i = 0; i < n; i++)				//LOOP TROUGH ALL THE VERTICES OF THE GRAPH
		{
			  //cout<<endl<<"STEP "<<i<<":\n________ ";
			  Minimum = INF;
		      for(int j = 0; j < n; j++){
                  if (T[j]==0){
                      if( Minimum > temp[j])
                      {
                          Minimum = temp[j];			//finding minimum L[s]
						  u = j;
                      }
			      }
			  }
			  //temp[u] = INFINITY;	//Assigning INIFINITY to the data structure allready visited to find the next minimum L
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
 */
};// end class

int main(int argc, char ** argv)
{
    string input, output;
    string outfile_all; // S.J. 03/30/2018 - to contain the output file for all subgraphs
    //cout<<argv[0]<<endl;
    parse_command_line(argc,argv,input,output,outfile_all); // S.J. - 03/30/2018 - added the function to parse the commanline arguments
    /*if(strcmp(argv[1],"-input")==0)
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
    }*/
    infile.open(input);
	outfile.open(output);
    
    //infile.open("/Users/cs4367/Documents/LouisLab/4OJI_adj.txt");
	//outfile.open("/Users/cs4367/Documents/LouisLab/4OJI_adj_out.txt");
    
    if(nods == 0 || nods == 1){ // 04/19/2018 - just print the information - this will be all 0s as nothing to be done for these
        
        outfile << "----------- Total number of blocks: " << numreg_blocks + numpseudo_blocks << endl;
        outfile << "----------- number of PK blocks: " << numpseudo_blocks << endl;
        outfile << "----------- number of regular blocks : " << numreg_blocks << endl;
        outfile << "-------------------------------------------------------------------------------------" << endl;
        
        infile.close();
        outfile.close();
        
        if(all_subgraphs){
            
            outfile.open(outfile_all);
            outfile << "----------- Total number of subgraphs between 2-9 vertices: 0" << endl;
            outfile << "----------- number of PK containing subgraphs: 0"  << endl;
            outfile << "----------- number of regular subgraphs : 0" << endl;
            outfile << "-------------------------------------------------------------------------------------" << endl;

            outfile.close();
        }
        
        return 0;
    }
    
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
        //int b_counter = 0; // S.J. 03/27/2018 - no need for this counter
        int Nmbrcmpnts=0; //we initialize the counter for the number of components
    
        if  (outfile.is_open()){
            //outfile << endl << "--------------------- Motif :" << motif_number << " -----------------------------" << endl;
        }
        else{
            cout<<"Cannot open outputfile!"<<endl;
            exit(0);
        }
        //outfile << " +++++++++++++++ number of vertices: " << nods << "    ====================" << endl;
        
        
        if(all_subgraphs){ // 04/03/2018 - to store articulation points
            
            artPoints.resize(nods);
        }
        
        
        Graph G(nods,edgs,dim,term);
        //while ( !infile.eof () && b_counter < 1){ // S.J. 03/27/2018 - no need for the while loop, statements executed only once
            
        G.init();
        //cout << endl << "ok";
        //b_counter ++; // S.J. 03/27/2018 - no need for the counter
		
        G.create ();
        //cout << endl << "ok2";
        //G.display();
        //cout << endl << "ok3";
		
        G.Bi_Connect(-1, 0);
        //outfile << endl << "----------- Summary information for Motif :" << motif_number <<"--------------------------------" << endl;
        outfile << "----------- Total number of blocks: " << numreg_blocks + numpseudo_blocks << endl;
        outfile << "----------- number of PK blocks: " << numpseudo_blocks << endl;
        outfile << "----------- number of regular blocks : " << numreg_blocks << endl;
        outfile << "-------------------------------------------------------------------------------------" << endl;
        //}
	
        infile.close();
        outfile.close();
        //system("pause");
        
        // S.J. 03/30/2018 - calculate all possible subgraphs and add it to the end of the list
        if(all_subgraphs){
            setup_artPoints(); // 04/03/2018 - to setup articulation points and its blocks
            //calc_possible_subgraphs();
            for(int i=0;i<Blocks.size();i++) // 04/03/2018 - calc all possible subgraphs starting from each block
                calc_possible_subgraphs(Blocks[i]);
            //cout << "All graph calculations done" << endl;
            print_subgraphs(outfile_all);
        }
    }
    else{
        cout<<"Cannot open inputfile!"<<endl;
        exit(0);
    }
    return 0;
}

/*int main()

{  

	infile.open("C:/Users/lpetn_000/desktop/RnaDv/2NOQ_B_matrix.txt");
    //infile.open("C:/users//desktop/RnaDv/swatimatrix.txt");
	//outfile.open("C:/users/petingi/desktop/RnaDv/.txt");
	outfile.open("C:/Users/lpetn_000/desktop/RnaDv/2NOQ_B_out.txt");

	clock_t start;
	clock_t end;
	float facto_time;
	float facto_irr_time;
	float facto_irr_e_time;
	string motif_number;
	cout<< endl << "Please enter the number of vertices: ";
	cin >>nods;
	
	cout<< endl << "Please enter motif number: ";
	cin >> motif_number;

	edgs=0;
	dim=0;
	term = 0;
    Visited=new bool [nods];
	int b_counter = 0;
	
    int Nmbrcmpnts=0; //we initialize the counter for the number of components
	
	outfile << endl << "--------------------- Motif :" << motif_number << " -----------------------------" << endl;


    //outfile << " +++++++++++++++ number of vertices: " << nods << "    ====================" << endl;
	Graph G(nods,edgs,dim,term);
	while ( !infile.eof () && b_counter < 1)
	{G.init();
	    //cout << endl << "ok";
	    b_counter ++;
		
	    G.create ();
		//G.display();
		 //cout << endl << "ok2";
		//G.display();
		//cout << endl << "ok3";
		


        G.Bi_Connect(-1, 0);
		outfile << endl << "----------- Summary information for Motif :" << motif_number << " --------------------------------" << endl;
		outfile << "----------- Total number of blocks: " << numreg_blocks + numpseudo_blocks << endl;
		outfile << "----------- number of PK blocks: " << numpseudo_blocks << endl;
		outfile << "----------- number of regular blocks : " << numreg_blocks << endl;
		outfile << "-------------------------------------------------------------------------------------" << endl;

	} 
	
    infile.close();
	outfile.close();
    system("pause");
  return 0;
 
}*/





