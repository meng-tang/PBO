#ifndef _EntropyMINIMIZE_H__
#define _EntropyMINIMIZE_H__
#include <BCD.h>
#include "PPBCEntropy.h"

void EntropyMinimize(int argc, char * argv[]);

// for Entropy-base image segmentation 
void EntropyMinimize(int argc, char * argv[])
{
	char * rootdir = "D:/repositories/PBO_public_code/";
	// parse arguments and set energy parameters
	char * imgname = argv[1];
	int numbin = atoi(argv[2]);
	double w_smooth = atof(argv[3]); // weight of smoothness term

	// load original image
	Image image(to_Cstr(rootdir<<"images/"<<imgname<<".bmp"),
		imgname,256/numbin,8); // 8 connect smoothness term
	image.print();
	int img_w = image.img_w;
	int img_h = image.img_h;
	srand( (unsigned int) time( NULL ) ); // RANDOM NUMBER INITIALIZER
	clock_t start,finish; // Timing

	// load bounding box
	Table2D<int> boximg = loadImage<RGB>(to_Cstr(rootdir<<"images/"<<imgname<<"_box.bmp"));
	int boxsize = countintable(boximg,0);
	outv(boxsize);
	Table2D<Label> initlabeling(img_w,img_h,NONE);
	Table2D<Label> hardconstraints(img_w,img_h,NONE);
	for(int x=0;x<img_w;x++){
		for(int y=0;y<img_h;y++){
			if(boximg[x][y]==0)//inside the box
				initlabeling[x][y] = OBJ;
			else{//outside the box
				initlabeling[x][y] = BKG;
				hardconstraints[x][y] = BKG;
			}
		}
	}
	Table2D<int> gt = loadImage<RGB>(to_Cstr(rootdir<<"images/"<<imgname<<"_groundtruth.bmp")); // ground truth

	cout<<"optimize with Block-coordinate-descent (GrabCut)"<<endl;
	start = clock();
	BCD bcd(image,w_smooth);
	bcd.initlabeling(initlabeling);
	bcd.hardconstraints = hardconstraints;
	bcd.optimize();
	finish = clock();
	double grabcuttime = (double)(finish-start)/CLOCKS_PER_SEC;

	cout<<"optimize with pPBC"<<endl;
	start = clock();
	PPBCEntropy ppbc = PPBCEntropy(image,w_smooth,1.0);
	ppbc.setpara(-4,4,0.1,20,true);
	Table2D<Label> ppbclabeling = initlabeling;
	double ppbcenergy = INFTY;
	double ppbctime;
	double ppbcerror;
	ppbc.hardconstraints = hardconstraints;
	int itr_count=0;
	start = clock();
	while(++itr_count)
	{
		cout<<itr_count<<" th iteration"<<endl;
		ppbc.setinitlabeling(ppbclabeling);
		ppbc.explore();
		BreakPoint best_bp = ppbc.SelectBestBP();
		cout<<"best breakpoint parameter: "<<best_bp.para<<endl;
		cout<<"best breakpoint energy: "<<best_bp.original_e<<endl;

		if(best_bp.original_e<ppbcenergy-0.1)
		{
			ppbclabeling = best_bp.solution;
			ppbcenergy = best_bp.original_e;
		}
		else
			break;
	}
	finish = clock();
	ppbctime = (double)(finish-start)/CLOCKS_PER_SEC;
	savebinarylabeling(image.img, ppbclabeling,
		to_Cstr(rootdir<<"images/"<<imgname<<"_PBO_withbox.bmp"));
	ppbcerror = geterrorrate(ppbclabeling, gt, boxsize);

	cout<<"optimize with pPBC (unsupervised)"<<endl;
	start = clock();
	ppbclabeling = Table2D<Label>(img_w,img_h,OBJ);
	double ppbcenergy_u = INFTY; // unsupervised (without bounding box)
	double ppbctime_u;
	double ppbcerror_u;
	itr_count=0;
	start = clock();
	ppbc.setpara(-4,10,0.1,20,true);
	while(++itr_count)
	{
		cout<<itr_count<<" th iteration"<<endl;
		ppbc.setinitlabeling(ppbclabeling);
		ppbc.explore();
		BreakPoint best_bp = ppbc.SelectBestBP();
		cout<<"best breakpoint parameter: "<<best_bp.para<<endl;
		cout<<"best breakpoint energy: "<<best_bp.original_e<<endl;

		if(best_bp.original_e<ppbcenergy_u-0.1)
		{
			ppbclabeling = best_bp.solution;
			ppbcenergy_u = best_bp.original_e;
		}
		else
			break;
	}
	finish = clock();
	ppbctime_u = (double)(finish-start)/CLOCKS_PER_SEC;
	ppbcerror_u = geterrorrate(ppbclabeling, gt, boxsize);

	// save resulting segmentation
	savebinarylabeling(image.img, bcd.current_labeling,
		to_Cstr(rootdir<<"images/"<<imgname<<"_BCD_withbox.bmp"));
	savebinarylabeling(image.img, ppbclabeling,
		to_Cstr(rootdir<<"images/"<<imgname<<"_PBO_withoutbox.bmp"));

	// energies
	cout<<"BCD energy: "<<bcd.current_e<<endl;
	cout<<"PBO energy: "<<ppbcenergy<<endl;
	cout<<"PBO energy (unsupervised): "<<ppbcenergy_u<<endl;

	// error rate
	double errorrate = geterrorrate(bcd.current_labeling, gt, boxsize);
	cout<<"BCD error rate: %"<<errorrate<<"%"<<endl;
	cout<<"PBO error rate: %"<<ppbcerror<<"%"<<endl;
	cout<<"PBO (unsupervised) error rate: %"<<ppbcerror_u<<"%"<<endl;

	// timing
	cout<<"BCD takes "<<grabcuttime<<" seconds!"<<endl;
	cout<<"PBO takes "<<ppbctime<<" seconds!"<<endl;
	cout<<"PBO (unsupervised) takes "<<ppbctime_u<<" seconds!"<<endl;
	return;
}

#endif