#pragma once
#include "Image.h"
#include <ezi/Cstr.h>
#include <entropy.h>

class BCD{
public:
	BCD(Image image_, double lambda_, double w_bits_ = 1, int max_iteration_ = 200, bool saveoption_ = false);
	void initlabeling(const Table2D<Label> & initiallabeling_);
	void updatemodel();
	bool updatelabeling();
	bool optimize(int max_iteration_ = 0);
	double trivialsolutionenergy();

	int max_iteration;
	Image image;
	int img_w;
	int img_h;

	int OBJsize;
	int BKGsize;
	vector<int> OBJhist;
	vector<int> BKGhist;

	Table2D<Label> hardconstraints;

	Table2D<bool> ROI;
	double init_e; // energy of initial labeling
	double current_e;
	Table2D<Label> current_labeling;
	// energy parameter
	double w_bits;
	double lambda; // weight of smoothness term.
	bool saveoption;
};

BCD::BCD(Image image_, double lambda_, double w_bits_, int max_iteration_, bool saveoption_)
	:image(image_),lambda(lambda_),max_iteration(max_iteration_),saveoption(saveoption_)
{
	img_w = image_.img_w;
	img_h = image_.img_h;
	w_bits = w_bits_;
	ROI = Table2D<bool>(img_w,img_h,true);
	hardconstraints = Table2D<Label>(img_w,img_h,NONE);
}

void BCD::initlabeling(const Table2D<Label> & initiallabeling_)
{
	current_labeling = initiallabeling_;
	ROI.reset(true);
	for(int y=0; y<img_h; y++)
	{
		for(int x=0; x<img_w; x++) 
		{ 
			if(current_labeling[x][y]==NONE)
				ROI[x][y] = false;
		}
	}
	init_e = getgrabcutenergy(image, w_bits, lambda, current_labeling);
	//cout<<"init energy: "<<init_e<<endl;
	current_e = init_e;
	updatemodel();
}
void BCD::updatemodel()
{
	OBJhist = vector<int>(image.colorbinnum,0);
	BKGhist = vector<int>(image.colorbinnum,0);
	OBJsize = 0;
	BKGsize = 0;
	for(int x =0;x<img_w;x++)
	{
		for(int y=0;y<img_h;y++)
		{
			if(current_labeling[x][y]==OBJ)
			{
				OBJsize++;
				OBJhist[image.colorlabel[x][y]]++;
			}
			else if(current_labeling[x][y]==BKG)
			{
				BKGsize++;
				BKGhist[image.colorlabel[x][y]]++;
			}
		}
	}
	//outv(OBJsize);
	//outv(BKGsize);
}

bool BCD::updatelabeling()
{
	if((OBJsize==0)||(BKGsize==0))
		return false;
	GraphType * g = new GraphType(img_w*img_h, 4*img_w*img_h); 
	g->add_node(img_w*img_h);    // adding nodes
	// add smoothness term
	addsmoothnessterm(g, image, lambda, ROI);
	// add unary term
	double source = 0,sink =0;
	for(int j=0;j<img_h;j++)
	{
		for(int i=0;i<img_w;i++)
		{
			int n = i+j*img_w;
			if(hardconstraints[i][j]==BKG) // hard constraints to background
				g->add_tweights(n,0,INFTY);
			else if(hardconstraints[i][j]==OBJ) // hard constraints to foreground
				g->add_tweights(n,INFTY,0);
			else if(ROI[i][j])
			{
				int idx = image.colorlabel[i][j];
				source  = max((double)BKGhist[idx],1.0)/(double)BKGsize;
				sink =max((double)OBJhist[idx],1.0)/(double)OBJsize;
				g->add_tweights(n,-log(source)/log(2.0)*w_bits,- log(sink)/log(2.0)*w_bits);//histogram natural log
			}
		}
	}
	double newflow = g -> maxflow();
	// set next labeling
	Table2D<Label> next_labeling(img_w,img_h,NONE);

	for (int y=0; y<img_h; y++) 
	{
		for (int x=0; x<img_w; x++) 
		{ 
			if(ROI[x][y])
			{
				int n = x+y*img_w;
				if(g->what_segment(n) == GraphType::SOURCE)
				{
					next_labeling[x][y]=OBJ;
				}
				else if(g->what_segment(n) == GraphType::SINK)
				{
					next_labeling[x][y]=BKG;
				}
			}
		}
	}
	delete g;
	g = NULL;
	double next_e = getgrabcutenergy(image, w_bits, lambda, next_labeling);
	if(next_e<current_e-0.1)
	{
		current_e = next_e;
		current_labeling = next_labeling;
		return true;
	}
	else
		return false;
}
bool BCD::optimize(int max_iteration_)
{
	int itr_count = 0;
	bool returnv = false;
	if(max_iteration_ !=0)
		max_iteration = max_iteration_;
	while(updatelabeling()&&(itr_count<max_iteration))
	{
		itr_count++;
		//outv(current_e);
		updatemodel();
		returnv = true;
		// assert(OBJsize!=0&&BKGsize!=0,"trivial solution for BCD!");
	}
	return returnv;
}

double BCD::trivialsolutionenergy()
{
	Table2D<Label> labeling(img_w,img_h,NONE);
	for (int y=0; y<img_h; y++) 
	{
		for (int x=0; x<img_w; x++) 
		{ 
			if(ROI[x][y])
				labeling[x][y] = OBJ;
		}
	}
	return getgrabcutenergy(image, w_bits, lambda, labeling);
}