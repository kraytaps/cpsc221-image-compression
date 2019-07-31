
/**
 *
 * toqutree (pa3)
 * significant modification of a quadtree .
 * toqutree.cpp
 * This file will be used for grading.
 *
 */

#include "toqutree.h"

toqutree::Node::Node(pair<int,int> ctr, int dim, HSLAPixel a)
	:center(ctr),dimension(dim),avg(a),NW(NULL),NE(NULL),SE(NULL),SW(NULL)
	{}

toqutree::~toqutree(){
	clear(root);
}

toqutree::toqutree(const toqutree & other) {
	root = copy(other.root);
}

toqutree & toqutree::operator=(const toqutree & rhs){
	if (this != &rhs) {
		clear(root);
		root = copy(rhs.root);
	}
	return *this;
}

toqutree::toqutree(PNG & imIn, int k){

/* This constructor grabs the 2^k x 2^k sub-image centered */
/* in imIn and uses it to build a quadtree. It may assume  */
/* that imIn is large enough to contain an image of that size. */

/* your code here */
	if (k == 0) {
		root = buildTree(&imIn, k);
	}

	else {

		PNG * square = new PNG(pow(2,k), pow(2,k));
		pair<int, int> kCtr(imIn.width() / 2, imIn.height() / 2);
		int rad = pow(2, k-1);
		pair<int,int> ul(kCtr.first - rad, kCtr.second - rad);
		pair<int,int> lr(kCtr.first + rad - 1, kCtr.second + rad - 1);

		int i = 0;
		int j;
		// cout<<"constructor"<<endl;
		// cout<<endl;
		for (int x = ul.first; x <= lr.first; x++){
			j = 0;
			for (int y = ul.second; y <= lr.second; y++){
				// cout<<i<<" "<<j<<"    "<<x<<" "<<y<<endl;
				*square->getPixel(i,j) = *imIn.getPixel(x,y);
				j++;
			}
			i++;
		}
		root = buildTree(square, k);
		delete square;
		square = NULL;
	}
}


int toqutree::sizeHelper(Node* croot) {

	if (croot == NULL) {
		return 0;
	}
	else if (croot->NW == NULL && croot->NE == NULL && croot->SW == NULL && croot->SE == NULL) {
		return 1;
	}
	else {
		return 1 + sizeHelper(croot->NW) + sizeHelper(croot->NW) +
								sizeHelper(croot->NW) + sizeHelper(croot->NW);
	}
}

int toqutree::size() {
/* your code here */
	return sizeHelper(root);
}


PNG * toqutree::buildBigPNG(PNG * im) {
// cout<<"IN BUILDBIGPNG"<<endl;
// cout<<endl;
	PNG * myPNG = new PNG(2*im->width(), 2*im->height());
	int i = 0;
	int j;


	for (unsigned int x = (im->width())/2; x < im->width(); x++) {
		j = 0;
		for (unsigned int y = (im->height())/2; y < im->height(); y++) {
			*myPNG->getPixel(i,j) = *im->getPixel(x,y);
			j++;
		}
		i++;
	}
	for (unsigned int x = 0; x < im->width(); x++) {
		j = 0;
		for (unsigned int y = (im->height())/2; y < im->height(); y++) {
			*myPNG->getPixel(i,j) = *im->getPixel(x,y);
			j++;
		}
		i++;
	}
	for (unsigned int x = 0; x < (im->width())/2; x++) {
		j = 0;
		for (unsigned int y = (im->height())/2; y < im->height(); y++) {
			*myPNG->getPixel(i,j) = *im->getPixel(x,y);
			j++;
		}
		i++;
	}

	i = 0;

	for (unsigned int x = (im->width())/2; x < im->width(); x++) {
		j = (im->height())/2;
		for (unsigned int y = 0; y < im->height(); y++) {
			*myPNG->getPixel(i,j) = *im->getPixel(x,y);
			j++;
		}
		i++;
	}
	for (unsigned int x = 0; x < im->width(); x++) {
		j = (im->height())/2;
		for (unsigned int y = 0; y < im->height(); y++) {
			*myPNG->getPixel(i,j) = *im->getPixel(x,y);
			j++;
		}
		i++;
	}
	for (unsigned int x = 0; x < (im->width())/2; x++) {
		j = (im->height())/2;
		for (unsigned int y = 0; y < im->height(); y++) {
			*myPNG->getPixel(i,j) = *im->getPixel(x,y);
			j++;
		}
		i++;
	}

	i = 0;

	for (unsigned int x = (im->width())/2; x < im->width(); x++) {
		j = (im->height())/2 + im->height();
		for (unsigned int y = 0; y < (im->height())/2; y++) {
			*myPNG->getPixel(i,j) = *im->getPixel(x,y);
			j++;
		}
		i++;
	}
	for (unsigned int x = 0; x < im->width(); x++) {
		j = (im->height())/2 + im->height();
		for (unsigned int y = 0; y < (im->height())/2; y++) {
			*myPNG->getPixel(i,j) = *im->getPixel(x,y);
			j++;
		}
		i++;
	}
	for (unsigned int x = 0; x < im->width()/2; x++) {
		j = im->height()/2 + im->height();
		for (unsigned int y = 0; y < im->height()/2; y++) {
			*myPNG->getPixel(i,j) = *im->getPixel(x,y);
			j++;
		}
		i++;
	}
	return myPNG;
}

//helper
double toqutree::getEntropy(PNG * bigim, PNG * im, stats * stat_im, int x, int y) {
	// cout<<"in getEntropy with k = "<<k<<endl;
	int center_sidelength = im->height()/2;
	pair<int, int> SE_ul(x, y);
	pair<int, int> SE_lr(x+center_sidelength-1, y+center_sidelength-1);
	pair<int, int> NW_lr(SE_ul.first - 1, SE_ul.second-1);
	pair<int, int> NW_ul(NW_lr.first - (center_sidelength-1), NW_lr.second - (center_sidelength-1));
	pair<int, int> NE_ul(SE_ul.first, NW_ul.second);
	pair<int, int> NE_lr(SE_lr.first, NW_lr.second);
	pair<int, int> SW_ul(NW_ul.first, SE_ul.second);
	pair<int, int> SW_lr(NW_lr.first, SE_lr.second);
	//
	// cout<<stat_im->entropy(SE_ul, SE_lr) << " " << stat_im->entropy(NE_ul, NE_lr) << " "
	// << stat_im->entropy(SW_ul, SW_lr) << " "<< stat_im->entropy(NW_ul, NW_lr)<<endl;
	return (stat_im->entropy(SE_ul, SE_lr) +
					stat_im->entropy(NE_ul, NE_lr) +
					stat_im->entropy(SW_ul, SW_lr) +
					stat_im->entropy(NW_ul, NW_lr)) / 4;
}


toqutree::Node * toqutree::buildTree(PNG * im, int k) {
/* your code here */

// Note that you will want to practice careful memory use
// In this function. We pass the dynamically allocated image
// via pointer so that it may be released after it is used .
// similarly, at each level of the tree you will want to
// declare a dynamically allocated stats object, and free it
// once you've used it to choose a split point, and calculate
// an average.

	stats * imgStats = new stats(*im);
	pair<int,int> ul(0,0);
	pair<int,int> lr(im->width()-1,im->height()-1);

	if (k == 0) {
		// cout<<"reached base case"<<endl;
		pair<int,int> nullctr(0, 0);
		Node * newNode = new Node(nullctr, k, imgStats->getAvg(ul, lr));


		delete imgStats;
		// cout<<"returning newNode k==0"<<endl;
		return newNode;

	}
	if (k == 1) {
		// cout<<"k==1"<<endl;
		pair<int,int> nullctr(1, 1);
		Node * newNode = new Node(nullctr, k, imgStats->getAvg(ul, lr));

		PNG * NWim = new PNG(1, 1);
		*NWim->getPixel(0,0) = *im->getPixel(0,0);
		PNG * NEim = new PNG(1, 1);
		*NEim->getPixel(0,0) = *im->getPixel(1,0);
		PNG * SWim = new PNG(1, 1);
		*SWim->getPixel(0,0) = *im->getPixel(0,1);
		PNG * SEim = new PNG(1, 1);
		*SEim->getPixel(0,0) = *im->getPixel(1,1);

		newNode->NW = buildTree(NWim, k-1);
		newNode->NE = buildTree(NEim, k-1);
		newNode->SW = buildTree(SWim, k-1);
		newNode->SE = buildTree(SEim, k-1);



		delete NWim;
		delete NEim;
		delete SWim;
		delete SEim;

		delete imgStats;
		return newNode;
	}

	else {
		PNG * bigPNG = buildBigPNG(im);

		PNG out(*im);
		if (k == 8) {
			out.writeToFile("images/myk8Image.png");
		}


		stats * bigStats = new stats(*bigPNG);
		int center_sidelength = im->width()/2;
		pair<int, int> ctr(0,0);
		double minEntropy = 1000000000.;

		for (unsigned int x = center_sidelength+center_sidelength/2; x < im->width() + center_sidelength/2 ; x++) {
			for (unsigned int y = center_sidelength + center_sidelength/2 ; y < im->height() + center_sidelength/2 ; y++) {
				double tmp_en = getEntropy(bigPNG, im, bigStats, x, y);
				if (tmp_en < minEntropy) {
					minEntropy = tmp_en;
					ctr.first = x;
					ctr.second = y;
				}
				cout<<"ctr: ("<<x<<", "<<y<<")\tte = "<<tmp_en<<"\tme = "<<minEntropy<<endl;
			}
		}
		Node * newNode = new Node(pair<int,int>(ctr.first - center_sidelength,
		ctr.second - center_sidelength), k, imgStats->getAvg(ul,lr));
		// cout<<endl;
		// cout<<endl;
		//
		// cout<<imgStats->getAvg(ul,lr)<<endl;
		// cout<<endl;

		pair<int, int> SE_ul(ctr.first, ctr.second);
		pair<int, int> SE_lr(ctr.first+center_sidelength-1, ctr.second+center_sidelength-1);
		pair<int, int> NW_lr(SE_ul.first-1, SE_ul.second-1);
		pair<int, int> NW_ul(NW_lr.first - (center_sidelength-1), NW_lr.second - (center_sidelength-1));
		pair<int, int> NE_ul(SE_ul.first, NW_ul.second);
		pair<int, int> NE_lr(SE_lr.first, NW_lr.second);
		pair<int, int> SW_ul(NW_ul.first, SE_ul.second);
		pair<int, int> SW_lr(NW_lr.first, SE_lr.second);

		PNG * NWim = new PNG((im->width())/2, (im->height())/2);

		int i = 0;
		int j;
		for (int x = NW_ul.first; x <= NW_lr.first; x++) {
			j = 0;
			for (int y = NW_ul.second; y <= NW_lr.second; y++) {
				*NWim->getPixel(i,j) = *bigPNG->getPixel(x,y);
				// cout<<"NW   "<<i<<" "<<j<<"   "<<x<<" "<<y<<endl;
				j++;
			}
			i++;
		}

		PNG * NEim = new PNG((im->width())/2, (im->height())/2);
		// cout<<"in makechildPNG"<<endl;
		// cout<<ul.first<<ul.second<<lr.first<<lr.second<<"   "<<ctr_ul.first<<ctr_ul.second<<ctr_lr.first<<ctr_lr.second<<endl;

		i = 0;
		for (int x = NE_ul.first; x <= NE_lr.first; x++) {
			j = 0;
			for (int y = NE_ul.second; y <= NE_lr.second; y++) {
				*NEim->getPixel(i,j) = *bigPNG->getPixel(x,y);
				// cout<<"NE   "<<i<<" "<<j<<"   "<<x<<" "<<y<<endl;
				j++;
			}
			i++;
		}

		PNG * SWim = new PNG((im->width())/2, (im->height())/2);
		// cout<<"in makechildPNG"<<endl;
		// cout<<ul.first<<ul.second<<lr.first<<lr.second<<"   "<<ctr_ul.first<<ctr_ul.second<<ctr_lr.first<<ctr_lr.second<<endl;

		i = 0;
		for (int x = SW_ul.first; x <= SW_lr.first; x++) {
			j = 0;
			for (int y = SW_ul.second; y <= SW_lr.second; y++) {
				*SWim->getPixel(i,j) = *bigPNG->getPixel(x,y);
				// cout<<"SW   "<<i<<" "<<j<<"   "<<x<<" "<<y<<endl;
				j++;
			}
			i++;
		}

		PNG * SEim = new PNG((im->width())/2, (im->height())/2);
		// cout<<"in makechildPNG"<<endl;
		// cout<<ul.first<<ul.second<<lr.first<<lr.second<<"   "<<ctr_ul.first<<ctr_ul.second<<ctr_lr.first<<ctr_lr.second<<endl;

		i = 0;
		for (int x = SE_ul.first; x <= SE_lr.first; x++) {
			j = 0;
			for (int y = SE_ul.second; y <= SE_lr.second; y++) {
				*SEim->getPixel(i,j) = *bigPNG->getPixel(x,y);
				// cout<<"SE   "<<i<<" "<<j<<"   "<<x<<" "<<y<<endl;
				j++;
			}
			i++;
		}
		newNode->NW = buildTree(NWim, k-1);
		newNode->NE = buildTree(NEim, k-1);
		newNode->SW = buildTree(SWim, k-1);
		newNode->SE = buildTree(SEim, k-1);

		delete NWim;
		delete NEim;
		delete SWim;
		delete SEim;

		delete bigPNG;
		delete bigStats;
		delete imgStats;
		return newNode;
	}
}


HSLAPixel toqutree::renderHelper(Node* node, int x, int y){
	if (node->NW == NULL && node->NE == NULL && node->SE == NULL && node->SW == NULL){
		return node->avg;
	}
	int dim = pow(2, node->dimension - 1);
	pair<int, int> loc(x, y);
	pair<int, int> nwCoord = NW(loc, node);
	pair<int, int> neCoord = NE(loc, node);
	pair<int, int> swCoord = SW(loc, node);
	pair<int, int> seCoord = SE(loc, node);

	if (nwCoord.first >= 0 && (nwCoord.first < dim) && nwCoord.second >= 0 && nwCoord.second < dim){
		return renderHelper(node->NW, nwCoord.first, nwCoord.second);
	}
	else if (neCoord.first >= 0 && (neCoord.first < dim) && neCoord.second >= 0 && neCoord.second < dim){
		return renderHelper(node->NE, neCoord.first, neCoord.second);
	}
	else if (swCoord.first >= 0 && (swCoord.first < dim) && swCoord.second >= 0 && swCoord.second < dim){
		return renderHelper(node->SW, swCoord.first, swCoord.second);
	}
	else if (seCoord.first >= 0 && (seCoord.first) < dim && swCoord.second >= 0 && swCoord.second < dim){
		return renderHelper(node->SE, seCoord.first, seCoord.second);
	}

}

pair<int, int> toqutree::NW(pair<int, int> loc, Node* subRoot) {
	int k = subRoot->dimension;
	int ctrx = subRoot->center.first;
	int ctry = subRoot->center.second;
	int parent = pow(2, k);
	int child = pow(2, k-1);
	int retX = (loc.first - ctrx + child + parent) % parent;
	int retY = (loc.second - ctry + child + parent) % parent;
	return make_pair(retX, retY);
}

pair<int, int> toqutree::NE(pair<int, int> loc, Node* subRoot) {
	int k = subRoot->dimension;
	int ctrx = subRoot->center.first;
	int ctry = subRoot->center.second;
	int parent = pow(2, k);
	int child = pow(2, k-1);
	int retX = (loc.first - ctrx + parent) % parent;
	int retY = (loc.second - ctry + child + parent) % parent;
	return make_pair(retX, retY);
}

pair<int, int> toqutree::SW(pair<int, int> loc, Node* subRoot) {
	int k = subRoot->dimension;
	int ctrx = subRoot->center.first;
	int ctry = subRoot->center.second;
	int parent = pow(2, k);
	int child = pow(2, k-1);
	int retX = (loc.first - ctrx + child + parent) % parent;
	int retY = (loc.second - ctry + parent) % parent;
	return make_pair(retX, retY);
}

pair<int, int> toqutree::SE(pair<int, int> loc, Node* subRoot) {
	int k = subRoot->dimension;
	int ctrx = subRoot->center.first;
	int ctry = subRoot->center.second;
	int parent = pow(2, k);
	// int child = pow(2, k-1);
	int retX = (loc.first - ctrx + parent) % parent;
	int retY = (loc.second - ctry + parent) % parent;
	return make_pair(retX, retY);
}

PNG toqutree::render(){

// My algorithm for this problem included a helper function
// that was analogous to Find in a BST, but it navigated the
// quadtree, instead.

/* your code here */
	int dim = root->dimension;
	int size = pow(2,dim);
	PNG output(size, size);

	for (int i = 0; i < size; i++){

		for (int j = 0; j < size; j++){
			*output.getPixel(i,j) = renderHelper(root, i, j);
		}
	}
	return output;

}

bool toqutree::isWithinTol(Node* root, double tol, HSLAPixel avg){
	// cout<<"is within tol begin"<<endl;
	if (root == NULL) {
		return false;
	}

	if (root->SE == NULL && root->SW == NULL && root->NE == NULL && root->NW == NULL){
		// cout<<"in if statement iswithintol"<<endl;
		return abs(root->avg.dist(avg)) <= tol;
	}
	else {
		// cout<<"in else iswinthineotl"<<endl;
		return isWithinTol(root->SE, tol, avg) & isWithinTol(root->SW, tol, avg) & isWithinTol(root->NE, tol, avg) & isWithinTol(root->NW, tol, avg);
	}
}
void toqutree::pruneHelper(Node* root, double tol){

	cout<<"in prune helper"<<endl;
		if (root == NULL){
			// cout<<"in base case prune"<<endl;
			return;
		}


		else if (isWithinTol(root->SE, tol, root->avg) &&
		isWithinTol(root->SW, tol, root->avg) &&
		isWithinTol(root->NE, tol, root->avg) &&
		isWithinTol(root->NW, tol, root->avg)){
			// cout<<"in else if case prune"<<endl;
			clearSubTrees(root);
		}

		else {
			// cout << "before pruning" << endl;
			pruneHelper(root->SE, tol);
			pruneHelper(root->SW, tol);
			pruneHelper(root->NE, tol);
			pruneHelper(root->NW, tol);
			// cout << "finish pruning" << endl;
		}
}
void toqutree::clearSubTrees(Node * root){
	cout << "before clear" << endl;
	clear(root->SE);
	clear(root->SW);
	clear(root->NE);
	clear(root->NW);
	root->SE = NULL;
	root->SW = NULL;
	root->NE = NULL;
	root->NW = NULL;
	cout << "after clear" << endl;
}
void toqutree::prune(double tol){
	pruneHelper(root,tol);

}
void toqutree::clearHelper(Node * & curr) {
	if (curr->NW == NULL && curr->NE == NULL && curr->SW == NULL && curr->SE == NULL) {
		delete curr;
	}
	else {
		clearHelper(curr->NW);
		clearHelper(curr->NE);
		clearHelper(curr->SW);
		clearHelper(curr->SE);
		delete curr;
	}
}

/* called by destructor and assignment operator*/
void toqutree::clear(Node * & curr){
/* your code here */
	clearHelper(curr);

}

/* done */
/* called by assignment operator and copy constructor */
toqutree::Node * toqutree::copy(const Node * other) {

/* your code here */

	if (other == NULL) {
		return NULL;
	}

	Node * thisRoot = new Node(other->center, other->dimension, other->avg);


	thisRoot->NW = copy(other->NW);
	thisRoot->NE = copy(other->NE);
	thisRoot->SW = copy(other->SW);
	thisRoot->SE = copy(other->SE);
	return thisRoot;
}
