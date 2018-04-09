#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <sstream>

using namespace std;

typedef vector<double> pnt;
typedef vector<vector<double>> points;
typedef vector<vector<vector<double>>> list_points;

struct intNode
{
	double data_median;
	bool is_leaf;
	pnt data_full_point; 		// level is depth%dimensions
	int level;
	intNode * l_intNode;
	intNode * r_intNode;
	points MBR;					//two points to define the rectangle(MBR)
	
};


//dimensions start from 0
struct comparator {
	int dim1;
    comparator(int dim) { this->dim1 = dim; }
    bool operator () (std::vector<double> i, std::vector<double> j) { return i[dim1]<j[dim1]; }
};

points sortPoints(points list, int dimension)
{	
	sort(list.begin(), list.end(), comparator(dimension));
	// cout<<"sort points"<<endl;
	return list;
}

void printPoints(points p)
{
	for(int i =0; i< p.size();i++){
		for (int j = 0; j < p[i].size(); j++)
		{
			cout<<p[i][j]<<' ';			
		}
		cout<<endl;
	}
}

class kdTree
{
public:
	intNode* root;
	int dimension;
	kdTree(int dim)
	{
		dimension = dim;
	}

	void buildStart(points &all_data)
	{
		//sort on all dimensions, get d lists of points
		// cout<<"build start"<<endl;
		list_points sorted_lists;
		for(int i=0;i<dimension;i++)
			sorted_lists.push_back(sortPoints(all_data, i));
		intNode* head = new intNode;
		// cout<<"end of buildStart1"<<endl;
		
		points rect;
		// cout<<"end of buildStart2"<<endl;

		for (int i = 0; i < 2; ++i)
		{
			pnt new_pnt;
			for (int j = 0; j < dimension; ++j)
			{
				if (i==0)				
					new_pnt.push_back(-1 * std::numeric_limits<float>::infinity());
				else	
					new_pnt.push_back(std::numeric_limits<float>::infinity());				
			}
			rect.push_back(new_pnt);
		}
		// cout<<"end of buildStart3"<<endl;

		head->MBR = rect;
		// cout<<"end of buildStart"<<endl;
		buildTree(sorted_lists, head, 0, rect);
		root = head;
		return;
	}

// Build the KD-tree
	void buildTree(list_points &sorted_lists, intNode* head, int level, points &rect)
	{
		//find median of levelth list -> access the size/2th element -> levelth 
		if (sorted_lists[level].size()==1)
		{
			//cout<<"size is 1 "<<endl;
			head->data_full_point = sorted_lists[0][0];
			head->is_leaf=true;
			head->data_median = sorted_lists[0][0][level];
			head->l_intNode=NULL;
			head->r_intNode=NULL;
			head->level=level;
			points rect;
			rect.push_back(head->data_full_point);
			rect.push_back(head->data_full_point);
			head->MBR = rect;
			cout<<head->data_median<<endl;
			printPoints(head->MBR);
			cout<<endl;
			// cout<<"leaf end"<<endl;
			return;
		}
		
		int no_of_points = sorted_lists[level].size();
		// cout<<"tot point " <<no_of_points<<" level : "<<level<<endl;
		double median = sorted_lists[level][floor((no_of_points-1)/2)][level];
		//cout<<"new node"<<endl;
		//cout<<median<<endl;
		head->data_median = median;
		head->is_leaf=false;
		head->MBR = rect;
		head->level = level;
		cout<<median<<endl;
		printPoints(head->MBR);
		cout<<endl;

		//cout<<"!"<<endl;
		
		// partition all other lists on basis of earlier median into two new lists 
		list_points left_lists;
		list_points right_lists;

		for(int i=0;i<sorted_lists.size();i++)
		{
			points left_points;			
			points right_points;
			// cout<<"1"<<endl;
			for (int j=0;j<no_of_points;j++)
			{
				// cout<<sorted_lists[i][j][level]<<endl;
				if(sorted_lists[i][j][level]<=median)
					{
						// cout<<"no of points"<<no_of_points<<endl;
						left_points.push_back(sorted_lists[i][j]);
					}
				else right_points.push_back(sorted_lists[i][j]);
				// cout<<"2"<<endl;
			}

			left_lists.push_back(left_points);
			right_lists.push_back(right_points);
		}

		
		// cout<<"left size is : "<<left_lists[0].size()<<endl;
		intNode* leftNode = new intNode;
		intNode* rightNode = new intNode;
		
		points left_rect;
		points right_rect;
		left_rect = rect;
		left_rect[1][head->level] = head->data_median;
		right_rect = rect;
		right_rect[0][head->level] = head->data_median;			

		buildTree(left_lists, leftNode, (level+1)%dimension , left_rect);
		head->l_intNode=leftNode;
		// cout<<"right size is : "<<right_lists[0].size()<<endl;
		buildTree(right_lists, rightNode, (level+1)%dimension, right_rect);
		head->r_intNode=rightNode;
		return;

	}
	~kdTree();
};

// // constructor and deconstructor have to be defined outside the class - o/w throws error in MacOS
// kdTree::kdTree()
// {
// 	//cout << "constructor" << endl;
// }

kdTree::~kdTree()
{
	//cout << "deconstructor" << endl;	
}


// min L2 distance of query point from MBR
/*
double distance_from_mbr(pnt data_point, points MBR)
{
	int dim = data_point.size();
	pnt delta;
	double dist = 0;
	for (int i = 0; i < dim; ++i)
	{
		if (data_point[i] < MBR[0][i])
		{
			delta.push_back(MBR[0][i] - data_point[i]);
		}
		else if(data_point[i] > MBR[1][i])
		{
			delta.push_back(data_point[i] - MBR[1][i]);
		}
		else delta.push_back(0);
	}

	for (int i = 0; i < dim; ++i)
	{
		dist = dist + delta[i]*delta[i]	;
	}
	dist = sqrt(dist);
	return dist;
}


// L2 distance of query point from another point
double distance_from_point(pnt query, pnt data_point)
{
	int dim = data_point.size();
	double dist = 0;
	for (int i = 0; i < dim; ++i)
	{
		double delta = abs(query[i] - data_point[i]);
		dist = dist + delta*delta;
	}
	dist = sqrt(dist);
	return dist;
}
*/


// // KNN query - best first
// void kNN_bestfirst(query_point, curr_node, min_dist)
// {


// }

// Read data points from dataset.txt
points readData(string dataset_file)
{
	// ifstream dataFile("dataset.txt");
	ifstream dataFile(dataset_file);
	string line;
	int dimension = 0;
	int numberOfPoints = 0;
	points all_points;


	if (dataFile.is_open())
    {
    	int i = 0;

    	while (! dataFile.eof() )
    	{
    		stringstream ss;
      		getline (dataFile,line);
      		ss.str(line);

      		if(line == ""){
      			continue;
      		}

    		if (i==0){
    			ss>>dimension;
    			ss>>numberOfPoints;
    			i++;
    			continue;
    		}     		
    		vector<double> point;
    		for(int j=0;j<dimension;j++){
    			double x;
    			ss>>x;
    			point.push_back(x);
    		}
    		all_points.push_back(point);
     		i++;
    	}
    }

    return all_points;
}


int main(int argc, char* argv[]) {

	// char* dataset_file = argv[1];
	string dataset_file = "dataset.txt";
	points all_points = readData(dataset_file);

	// [TODO] Construct kdTree using dataset_file here	
	kdTree mykdtree(all_points[0].size());
	mykdtree.buildStart(all_points);
	// mykdtree.set_MBRs(kdroot);

	// Request name/path of query_file from parent by just sending "0" on stdout
	// cout << 0 << endl;

	// Wait till the parent responds with name/path of query_file and k | Timer will start now
	// char* query_file = new char[100];
	// int k;
	// cin >> query_file >> k;
	// cerr << dataset_file << " " << query_file << " " << k << endl;

	// [TODO] Read the query point from query_file, do kNN using the kdTree and output the answer to results.txt

	// Convey to parent that results.txt is ready by sending "1" on stdout | Timer will stop now and this process will be killed
	cout << 1 << endl;
}