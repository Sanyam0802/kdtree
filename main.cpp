#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <math.h>

using namespace std;

typedef vector<vector<double>> points;
typedef vector<vector<vector<double>>> list_points;

struct leaf
{
	vector <double> data;
};

struct intNode
{
	double data;
	// level is depth%dimensions
	int level;
	intNode * l_intNode;
	intNode * r_intNode;
	leaf * l_leaf;
	leaf * r_leaf;
};


//dimensions start from 0
struct comparator {
    comparator(int dim) { this->dim = dim; }
    bool operator () (std::vector<double> i, std::vector<double> j) { return i[dim]<j[dim]; }
    int dim;
};

points sortPoints(points list, int dimension)
{
	sort(list.begin(), list.end(), comparator(dimension));
	return list;
}

class kdTree
{
public:
	intNode* root;
	int dimension;
	kdTree();

	void buildStart(points all_data)
	{
		//sort on all dimensions, get d lists of points
		list_points sorted_lists;
		for(int i=0;i<dimension;i++)
			sorted_lists.push_back(sortPoints(all_data, i));
		buildTree(sorted_lists, 0);
	}
	intNode* buildTree(list_points sorted_lists, int level)
	{
	//find median of levelth list -> access the size/2th element -> levelth 
		int no_of_points = sorted_lists[level].size();
		double median = sorted_lists[level][ceil(no_of_points/2)][level];
	//create intNode with this data
		intNode* newNode;
		newNode->data = median;
	// partition all other lists on basis of earlier median into two new lists 
		list_points left_lists;
		list_points right_lists;
		for(int i=0;i<sorted_lists.size();i++)
		{
			points left_points;
			points right_points;
			for (int j=0;j<no_of_points;j++)
			{
				if(sorted_lists[i][j][level]<=median)
					left_points.push_back(sorted_lists[i][j]);
				else right_points.push_back(sorted_lists[i][j]);
			}
			left_lists.push_back(left_points);
			right_lists.push_back(right_points);
		}
	//if size of lists is 1 after partitioning then 
		//make leaf nodes of this point as children, return this intNode with its children intNodes as null
	//else
		//call buildtree of left and right partition with level+1
		//return jo hua voh left aur right child ban gaya
		//make its leaf nodes null
		//return intNode jo banaya th

		if(left_lists[0].size()==1)
		{
			leaf* left_child;
			left_child->data = left_lists[0][0];
			newNode->l_leaf = left_child;
			newNode->l_intNode = NULL;
		} 
		else
		{
			newNode->l_leaf=NULL;
			newNode->l_intNode = buildTree(left_lists, level+1);
		}
		if (right_lists[0].size()==1)
		{
			leaf* right_child;
			right_child->data = right_lists[0][0];
			newNode->r_leaf = right_child;
			newNode->r_intNode = NULL;
		}
		else
		{
			newNode->r_leaf=NULL;
			newNode->r_intNode = buildTree(right_lists, level+1);
		}
		return newNode;

	}
	~ kdTree();
};





int main(int argc, char* argv[]) {

	char* dataset_file = argv[1];

	// [TODO] Construct kdTree using dataset_file here



	// Request name/path of query_file from parent by just sending "0" on stdout
	cout << 0 << endl;

	// Wait till the parent responds with name/path of query_file and k | Timer will start now
	char* query_file = new char[100];
	int k;
	cin >> query_file >> k;
	// cerr << dataset_file << " " << query_file << " " << k << endl;

	// [TODO] Read the query point from query_file, do kNN using the kdTree and output the answer to results.txt

	// Convey to parent that results.txt is ready by sending "1" on stdout | Timer will stop now and this process will be killed
	cout << 1 << endl;
}
