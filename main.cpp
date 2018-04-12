#include <iostream>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <math.h>
#include <fstream>
#include <sstream>
#include <queue>
#include <chrono>
#include <list>

using namespace std;

typedef vector<double> pnt;
typedef vector<vector<double>> points;
typedef vector<vector<vector<double>>> list_points;

struct intNode
{
	double data_median;
	double dist_to_query = std::numeric_limits<float>::infinity();
	bool is_leaf;
	pnt data_full_point; 		// level is depth%dimensions
	int level;
	intNode * l_intNode;
	intNode * r_intNode;
	points MBR;					//two points to define the rectangle(MBR)
	bool not_null = true;
	
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
		list_points sorted_lists;
		for(int i=0;i<dimension;i++)
			sorted_lists.push_back(sortPoints(all_data, i));
		intNode* head = new intNode;		
		points rect;
		// cout<<"sorted"<<endl;
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

		head->MBR = rect;
		buildTree(sorted_lists, head, 0, rect);
		// cout<<"hope it reaches here"<<endl;
		root = head;
		sorted_lists.clear();
		sorted_lists.shrink_to_fit();
		return;
	}

	// Build the KD-tree
	void buildTree(list_points &sorted_lists, intNode* head, int level, points &rect)
	{
		//find median of levelth list -> access the size/2th element -> levelth 
		// cout<<sorted_lists[level].size()<<endl;

		if (sorted_lists[level].size()==0)
		{
			head->not_null = false;
			head->is_leaf = false;
			// cout<<"eawda"<<endl;
			return;
		}

		if (sorted_lists[level].size()==1)
		{
			// cout<<"calledddd"<<endl;
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
			// cout<<head->data_median<<endl;
			// printPoints(head->MBR);
			// cout<<endl;
			return;
		}

		
		int no_of_points = sorted_lists[level].size();
		double median = sorted_lists[level][floor((no_of_points-1)/2)][level];
		head->data_median = median;
		head->is_leaf=false;
		head->MBR = rect;
		head->level = level;
		// cout<<"called"<<endl;

		// cout<<median<<endl;
		// printPoints(head->MBR);
		// cout<<endl;

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
					{
						left_points.push_back(sorted_lists[i][j]);
					}
				else right_points.push_back(sorted_lists[i][j]);
			}

			left_lists.push_back(left_points);
			right_lists.push_back(right_points);
		}

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

		left_lists.clear();
		left_lists.shrink_to_fit();

		buildTree(right_lists, rightNode, (level+1)%dimension, right_rect);
		head->r_intNode=rightNode;

		right_lists.clear();
		right_lists.shrink_to_fit();

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

double distance_from_mbr(pnt& data_point, points& MBR)
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
double distance_from_point(pnt& query, pnt& data_point)
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


struct comparator_max_heap {
	// pnt query;
 //    comparator_max_heap(pnt& query) { this->query = query; }
    bool operator () (pair<pnt, double> const& p1, pair<pnt, double> const& p2) 
    { 
    	return p1.second < p2.second; 
    }
};


struct comparator_min_heap {
    bool operator () (const intNode* n1, const intNode* n2) 
    { 
    	return n1->dist_to_query > n2->dist_to_query;
    }
};

// KNN query - best first
priority_queue< pair<pnt,double>, vector<pair<pnt,double>>, comparator_max_heap> kNN_bestfirst(int k, pnt& query_point, intNode* head, points& all_points)
{
	// answer set:	initialise Max Heap of pnts of size k with k random points
	// candidate:	initilise Min Heap with root MBR
	// while MBR is closer to query than top node in answer_set OR candidate is not empty
		//	pop top MBR
		// if MBR is a leaf AND closer to query than top node in answer set
			// pop top of answer set and insert node in answer set
		// else
			// insert those children of MBR which are closer to query point than the top node in the answer set

	priority_queue< pair<pnt,double>, vector<pair<pnt,double>>, comparator_max_heap> answer_set;
	priority_queue<intNode*, vector<intNode*>, comparator_min_heap> candidate;
	// initialise wih first k points
	for (int i = 0; i < k; ++i)
	{
		answer_set.push(make_pair(all_points[i], distance_from_point(query_point, all_points[i])));
	}
	
	
	head->dist_to_query = distance_from_mbr(query_point, head->MBR);

	candidate.push(head);
	// cout<<"WHILE . . . "<<endl;
	while(!candidate.empty() && (candidate.top()->dist_to_query < answer_set.top().second) )
	{
		// cout<<"candidate top median, level: "<<candidate.top()->data_median<< ", "<< candidate.top()->level<<endl;
		intNode* top_MBR = candidate.top();
		// cout<< "MBR. . . . "<<endl;
		// printPoints(top_MBR->MBR);
		candidate.pop();
		if (top_MBR->is_leaf && (answer_set.top().second > top_MBR->dist_to_query))
		{	
			// cout<<"IF. . "<<endl;
			answer_set.pop();
			answer_set.push(make_pair(top_MBR->data_full_point, distance_from_point(query_point, top_MBR->data_full_point)));
			// cout<<"size answers set: "<<answer_set.size()<<endl;
		}	
		else
		{
			// cout<<"ELSE. .  "<<endl;
			intNode* leftChild = top_MBR->l_intNode;
			intNode* rightChild = top_MBR->r_intNode;
			// cout<<"defined l, r child"<<endl;
			double right_distance;
			if (rightChild->not_null == false){
				// cout<<"yaya"<<endl;
				right_distance = std::numeric_limits<float>::infinity();
			}
			else
				right_distance = distance_from_mbr(query_point, rightChild->MBR);
			// cout<<". . "<<endl;
			double left_distance = distance_from_mbr(query_point, leftChild->MBR);
			// cout<<"calculated distance of left and right MBRs from query_point"<<endl;

			

			leftChild->dist_to_query = left_distance;
			rightChild->dist_to_query = right_distance;
			// cout<<"set dist_to_query"<<endl;

			// cout<<"left distance: "<<left_distance<<",  Right distance: "<<right_distance<<endl<<endl;

			if (left_distance < answer_set.top().second)
			{
				candidate.push(leftChild);
			}
			if (right_distance < answer_set.top().second)
			{
				candidate.push(rightChild);
			}
		}
	}
	return answer_set;
}

priority_queue< pair<pnt,double>, vector<pair<pnt,double>>, comparator_max_heap> kNN_sequential_scan(int k, pnt& query_point, intNode* head, points& all_points)
{
	priority_queue< pair<pnt,double>, vector<pair<pnt,double>>, comparator_max_heap> answer_set;
	for (int i = 0; i < k; ++i)
	{
		answer_set.push(make_pair(all_points[i], distance_from_point(query_point, all_points[i])));
	}
	for (int i = k; i < all_points.size(); ++i)
	{	
		double dist = distance_from_point(query_point, all_points[i]);
		if (answer_set.top().second > dist)
		{
			answer_set.pop();
			answer_set.push(make_pair(all_points[i], dist));
		}
	}

	return answer_set;

}


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

void write_result(priority_queue< pair<pnt,double>, vector<pair<pnt,double>>, comparator_max_heap> answer_set)
{	
	list<pnt> result;
	int dim = answer_set.top().first.size();
	while(!answer_set.empty())
	{
		result.push_back(answer_set.top().first);
		answer_set.pop();
	}

	ofstream results_file;
	results_file.open("results.txt");
	for (int i = 0; i < answer_set.size(); ++i)
	{
		vector<double> res = result.back();
		for (int j = 0; i < dim; ++i)
		{
			results_file<< res[j]<<" ";
		}
		result.pop_back();
		results_file<<"\n";
	}
	return;
}

int main(int argc, char* argv[]) {

	char* dataset_file = argv[1];
	char* queries_file = argv[2];
	char* K_neighbours = argv[3];
	int K = std::atoi(K_neighbours);

	// string dataset_file = "dataset.txt";
	points all_points = readData(dataset_file);
	// cout<<"reached"<<endl;
	// [TODO] Construct kdTree using dataset_file here	
	kdTree mykdtree(all_points[0].size());
	mykdtree.buildStart(all_points);
	cout<<"Build tree completed"<<endl;
	points all_query_points = readData(queries_file);
	cout<<"points loaded"<<endl;
	exit(1);
	priority_queue< pair<pnt,double>, vector<pair<pnt,double>>, comparator_max_heap> answer_set;

	auto start_seq = std::chrono::system_clock::now();
	for(int i = 0; i< all_query_points.size();i++)
	{
		// answer_set =  kNN_bestfirst(K, all_query_points[i], mykdtree.root, all_points);
		answer_set =  kNN_sequential_scan(K, all_query_points[i], mykdtree.root, all_points);
		// while(!answer_set.empty())
		// {
		// 	cout<<answer_set.top().second<<endl;
		// 	answer_set.pop();
		// }
		// cout<< "kNN run kiya humne yayy"<<endl<<endl;

	}
	auto end_seq = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds_seq = end_seq- start_seq;
	cout<<K<<" nearest neighbour query"<<endl;
	cout<<all_query_points[0].size()<<" dimension"<<endl;
	cout<<all_query_points.size()<<" query points"<<endl;
	cout<<all_points.size()<<" points in tree construction"<<endl<<endl;

	cout<< "Sequential";
	cout<< "Elapsed time: "<<elapsed_seconds_seq.count() << " s"<<endl;
	double avg_time_seq = elapsed_seconds_seq.count()/(1.0*all_query_points.size());
	cout<< "Average Elapsed time per query: "<<avg_time_seq << " s"<<endl;

	cout<<endl<<endl;

	auto start_best = std::chrono::system_clock::now();
	for(int i = 0; i< all_query_points.size();i++)
	{
		// for(int j = 0; j < all_query_points[i].size();j++){
		// 	cout<<all_query_points[i][j]<<' ';
		// }
		cout<<i<<endl;
		answer_set =  kNN_bestfirst(K, all_query_points[i], mykdtree.root, all_points);
		// answer_set =  kNN_sequential_scan(K, all_query_points[i], mykdtree.root, all_points);
		// while(!answer_set.empty())
		// {
		// 	cout<<answer_set.top().second<<endl;
		// 	answer_set.pop();
		// }
		// cout<< "kNN run kiya humne yayy"<<endl<<endl;

	}
	auto end_best = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds_best = end_best - start_best;
	cout<<"Best First";
	cout<< "Elapsed time: "<<elapsed_seconds_best.count() << "s"<<endl;
	double avg_time_best = elapsed_seconds_best.count()/(1.0*all_query_points.size());
	cout<< "Average Elapsed time per query: "<<avg_time_best << "s"<<endl;

	cout<<endl<<"Averge time in seq to best "<<avg_time_seq/avg_time_best<<endl;

	// for (int i = 0; i< all_points.size();i++){
	// 	cout<<distance_from_point(v1,all_points[i])<<endl;
	// }
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
