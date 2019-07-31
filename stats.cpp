
#include "stats.h"



stats::stats(PNG & im){
/* your code here */
  double mult = PI/180;
  sumHueX.resize(im.width(), vector<double>(im.height(), 0));
  sumHueY.resize(im.width(), vector<double>(im.height(), 0));
  sumSat.resize(im.width(), vector<double>(im.height(), 0));
  sumLum.resize(im.width(), vector<double>(im.height(),0));
  hist.resize(im.width(), vector<vector<int>> (im.height(), vector<int>(36,0)));
  for (unsigned int i = 0; i < (unsigned int) im.width(); i++){
    for (unsigned int j = 0; j < (unsigned int) im.height(); j++){
      HSLAPixel curr = *im.getPixel(i,j);
      if (i == 0 && j == 0){
        sumHueX[i][j] = cos(mult*curr.h);
        sumHueY[i][j] = sin(mult*curr.h);
        sumSat[i][j] = curr.s;
        sumLum[i][j] = curr.l;
        hist[i][j][curr.h/10] = 1;
      }
      else if (i != 0 && j == 0){
        sumHueX[i][j] = cos(mult*curr.h) + sumHueX[i-1][j];
        sumHueY[i][j] = sin(mult*curr.h) + sumHueY[i-1][j];
        sumSat[i][j] = curr.s + sumSat[i-1][j];
        sumLum[i][j] = curr.l + sumLum[i-1][j];
        hist[i][j] = hist[i-1][j];
        hist[i][j][curr.h/10]++;
      }
      else if (i == 0 && j != 0){
        sumHueX[i][j] = cos(mult*curr.h) + sumHueX[i][j-1];
        sumHueY[i][j] = sin(mult*curr.h) + sumHueY[i][j-1];
        sumSat[i][j] = curr.s + sumSat[i][j-1];
        sumLum[i][j] = curr.l + sumLum[i][j-1];
        hist[i][j] = hist[i][j-1];
        hist[i][j][curr.h/10]++;
      }
      else {
        sumHueX[i][j] = cos(mult*curr.h) + sumHueX[i][j-1] + sumHueX[i-1][j] - sumHueX[i-1][j-1];
        sumHueY[i][j] = sin(mult*curr.h) + sumHueY[i][j-1] + sumHueY[i-1][j] - sumHueY[i-1][j-1];
        sumSat[i][j] = curr.s + sumSat[i][j-1] + sumSat[i-1][j] - sumSat[i-1][j-1];
        sumLum[i][j] = curr.l + sumLum[i][j-1] + sumLum[i-1][j] - sumLum[i-1][j-1];
        for (int k = 0; k < 36; k++) {
          hist[i][j][k] = hist[i][j-1][k] + hist[i-1][j][k] - hist[i-1][j-1][k];
        }
        hist[i][j][curr.h/10]++;
      }
    }
  }
}

long stats::rectArea(pair<int,int> ul, pair<int,int> lr){
  /* your code here */
  return abs(((lr.first-ul.first)+1) * ((lr.second-ul.second)+1));
}

double stats::getPxArea(pair<int,int> ul, pair<int,int> lr, vector<vector<double>> * type){
  if (ul.first == 0 && ul.second == 0)
    return (*type)[lr.first][lr.second];
  else if (ul.first != 0 && ul.second == 0)
    return (*type)[lr.first][lr.second] - (*type)[ul.first-1][lr.second];
  else if (ul.first == 0 && ul.second != 0)
    return (*type)[lr.first][lr.second] - (*type)[lr.first][ul.second-1];
  else
    return (*type)[lr.first][lr.second] - (*type)[ul.first-1][lr.second] - (*type)[lr.first][ul.second-1] + (*type)[ul.first-1][ul.second-1];
}
HSLAPixel stats::getAvg(pair<int,int> ul, pair<int,int> lr){

/* your code here */
  double mult = 180/PI;
  HSLAPixel avgPix;
  if (getPxArea(ul, lr, &sumHueY) < 0)
    avgPix.h = atan2(getPxArea(ul,lr, &sumHueY)/rectArea(ul,lr), getPxArea(ul,lr,&sumHueX)/rectArea(ul,lr))*mult + 360;
  else
    avgPix.h = atan2(getPxArea(ul,lr, &sumHueY)/rectArea(ul,lr), getPxArea(ul,lr,&sumHueX)/rectArea(ul,lr))*mult;
  avgPix.s = getPxArea(ul, lr, &sumSat) / rectArea(ul, lr);
  avgPix.l = getPxArea(ul, lr, &sumLum) / rectArea(ul, lr);
  avgPix.a = 1.0;
  return avgPix;
}

vector<int> stats::buildHist(pair<int,int> ul, pair<int,int> lr){

    vector<int> recthist(36);

    if (ul.first == 0 && ul.second == 0) {
      for (int k = 0; k < 36; k++) {
        recthist[k] = hist[lr.first][lr.second][k];
      }
    }

    else if (ul.first == 0 && ul.second != 0) {
      for (int k = 0; k < 36; k++) {
        recthist[k] = hist[lr.first][lr.second][k] - hist[lr.first][ul.second-1][k];
      }
    }

    else if (ul.first != 0 && ul.second == 0) {
      for (int k = 0; k < 36; k++) {
        recthist[k] = hist[lr.first][lr.second][k] - hist[ul.first-1][lr.second][k];
      }
    }

    else {
      for (int k = 0; k < 36; k++) {
        recthist[k] = hist[lr.first][lr.second][k] - hist[lr.first][ul.second-1][k] -
                      hist[ul.first-1][lr.second][k] + hist[ul.first-1][ul.second-1][k];
      }
    }
    return recthist;
}

// takes a distribution and returns entropy
// partially implemented so as to avoid rounding issues.
double stats::entropy(vector<int> & distn,int area){

    double entropy = 0.;

    for (int i = 0; i < 36; i++) {
        if (distn[i] > 0 )
            entropy += ((double) distn[i]/(double) area)
                                    * log2((double) distn[i]/(double) area);
    }

    return  -1 * entropy;

}

double stats::entropy(pair<int,int> ul, pair<int,int> lr){

  vector<int> recthist = buildHist(ul, lr);
  int rectangle = (int) rectArea(ul,lr);
  double entropyVal = entropy(recthist, rectangle);
  return entropyVal;
  /* your code here */

}
