#ifndef processV3_H_
#define processV3_H_

#include "Bmath/pnorm.c"

using namespace std;


void PairWiseRRCal(In3str1v * paraFA04, Para_18 * para_18 , map <llong, vector <BaseType>  >  & SNPList ,  StarRsult  *  All_Stat[])
{
	double  CalResult;
	llong Dis=0;
	statementVar  Var;
	map<llong, vector <BaseType> >  :: iterator key2;
	map<llong, vector <BaseType> >  :: iterator key2_se;

	key2=SNPList.begin();
	Var.Asize= (key2->second).size();
	for ( ; key2!=SNPList.end(); key2++)
	{
		key2_se=key2 ; key2_se++ ;
		for( ; key2_se!=SNPList.end(); key2_se++)
		{
			Dis=(key2_se->first)-(key2->first);
			if ( Dis> (paraFA04->WindowSize))
			{
				break ;
			}
			if (cal_RR_MA( key2->second , key2_se->second  ,CalResult, Var) ==1)
			{
				int count=int((key2->first)/(paraFA04->bin));
				int SeCount=int(Dis/(paraFA04->bin));
				All_Stat[SeCount][count].Count++;
				All_Stat[SeCount][count].sumRR+=CalResult;
			}
		}
	}
	//	return 1;
}



void DrawTowDepthSliding( string File , In3str1v * paraFA04, Para_18 * para_18 ,string IDA, string IDB)
{

	//OUT<<"#Chr\tStart\tEnd\tMean_r^2_"<<IDA<<"\tSum_r^2_"<<IDA<<"\tCount_"<<IDA<<"\tMean_r^2_"<<IDB<<"\tSum_r^2_"<<IDB<<"\tCount_"<<IDB<<"\tMeanRRDiff("<<IDA<<"-"<<IDB<<")<<ZScore\tPvalue\n\n";


	string OutPlotTmp=(paraFA04->InStr2)+".tmp";
	ofstream  OUTPLOT (OutPlotTmp.c_str());

	igzstream LIST (File.c_str(),ifstream::in);

	string chr="";
	getline(LIST,chr);
	getline(LIST,chr);
	getline(LIST,chr);
	getline(LIST,chr);
	string Start,End,MeanRRV1,SumESV1,CountV1,MeanRRV2;

	map <string , int >  VecChr;
	map <string , map <int , string > >  ValueV1;
	map <string , map <int , string > >  ValueV2;
	map <string , map <int , string > > :: iterator iTValue;

	map  <string , int >  :: iterator itVecChr ;
	while(!LIST.eof())
	{
		string  line ;
		getline(LIST,line);
		if (line.length()<=0) { continue ;}
		istringstream isone (line,istringstream::in);
		isone>>chr>>Start>>End>>MeanRRV1>>SumESV1>>CountV1>>MeanRRV2 ;
		itVecChr=VecChr.find(chr);
		if (itVecChr==VecChr.end())
		{
			VecChr.insert(map <string,int>::value_type(chr,1));
			map <int , string >  DataEEV1;
			map <int , string >  DataEEV2;
			DataEEV1[0]=MeanRRV1;
			DataEEV2[0]=MeanRRV2;
			ValueV1.insert(map <string , map <int , string > > ::value_type(chr,DataEEV1));
			ValueV2.insert(map <string , map <int , string > > ::value_type(chr,DataEEV2));
		}
		else
		{
			iTValue=ValueV1.find(chr);
			(iTValue->second).insert(map <int , string >::value_type(itVecChr->second,MeanRRV1));
			iTValue=ValueV2.find(chr);
			(iTValue->second).insert(map <int , string >::value_type(itVecChr->second,MeanRRV2));
			itVecChr->second++;
		}
	}
	LIST.close();

	int MaxXX=0;
	OUTPLOT<<"Site";
	for(itVecChr=VecChr.begin();itVecChr!=VecChr.end();itVecChr++)
	{
		if (MaxXX<(itVecChr->second)) {MaxXX=itVecChr->second;}
		OUTPLOT<<"\t"<<IDA<<"_"<<itVecChr->first<<"\t"<<IDB<<"_"<<itVecChr->first;
	}
	OUTPLOT<<endl;


	map <int , string > ::  iterator  itEEEE;
	for (int ii=0 ; ii< MaxXX ; ii++)
	{
		OUTPLOT<<ii;
		for(itVecChr=VecChr.begin();itVecChr!=VecChr.end();itVecChr++)
		{
			iTValue=ValueV1.find(chr);
			itEEEE=(iTValue->second).find(ii);
			if (itEEEE==(iTValue->second).end())
			{
				OUTPLOT<<"\tNA";
			}
			else
			{
				OUTPLOT<<"\t"<<itEEEE->second ;
			}

			iTValue=ValueV2.find(chr);
			itEEEE=(iTValue->second).find(ii);
			if (itEEEE==(iTValue->second).end())
			{
				OUTPLOT<<"\tNA";
			}
			else
			{
				OUTPLOT<<"\t"<<itEEEE->second ;
			}
		}
		OUTPLOT<<endl;

	}


	OUTPLOT.close();

	string OutPlotr=(paraFA04->InStr2)+".r";
	ofstream  OUTR (OutPlotr.c_str());

	OUTR<<""
		"\n"
		"library(ggplot2);\nlibrary(reshape2);\n"
		"read.table(\""<<OutPlotTmp<<"\",header=T)->data;\n"
		"data2<-melt(data,id=1);\n"
		"p <- ggplot(data = data2)+geom_tile(aes(x = Site,y = variable, fill=value))+scale_fill_gradient(limit=c(0,1),low ='white',high ='red',na.value=\"grey50\")\n"
		"pdf(\""<<(paraFA04->InStr2)<<".pdf\",h=6,w=25);\n"
		"p\n"
		"dev.off()\n"
		"png(\""<<(paraFA04->InStr2)<<".png\",h=6,w=25);\n"
		"p\n"
		"dev.off()\n"<<endl;
	OUTR.close();

	char   buf[2048]={'\0'};
	string cc="which  Rscript  2> /dev/null ";
	memset( buf, '\0', sizeof(buf) );
	FILE   *stream ;
	stream=popen(cc.c_str(),"r") ;
	fread( buf, sizeof(char), sizeof(buf), stream);
	string binPath=buf;
	binPath=binPath.substr(0,binPath.length()-1);
	if (binPath == "" )
	{
		cout <<"\twarning: can't find the [Rscript] in your $PATH ; no png Figure Out"<<endl;
		cout <<"\t\tRscript "<<OutPlotr<<endl;
	}
	else
	{
		if (paraFA04->TF)
		{
			//cc=binPath+"\t"+OutPlotr ;
			cc=binPath+"\t"+OutPlotr+"  ; rm  -rf  "+ OutPlotr ;
			std::system(cc.c_str()) ;
		}
		else
		{
			cc=binPath+"\t"+OutPlotr ;
			std::system(cc.c_str()) ;
			cout <<"\t\tRePlot : Rscript  "<<OutPlotr<<endl;
		}
	}
}



int SilidingRRCalTowGroup( In3str1v * paraFA04, Para_18 * para_18 , map <string,map <llong, vector <BaseType> > >  & SNPListV1 ,   map <string,map <llong, vector <BaseType> > >  & SNPListV2 , map <string,llong>  & MaxSNPSite  ,   map <string,int >  & GroupID )
{
	string Stat=(paraFA04->InStr2)+".winRR.temp.gz";
	ogzstream OUT ((Stat).c_str());
	if((!OUT.good()))
	{
		cerr << "open OUT File error: "<<Stat<<endl;
		return 0;
	}

	string IDA="A";
	string IDB="B";
	map <string, int> :: iterator  ITID ;
	for (ITID=GroupID.begin(); ITID!=GroupID.end();ITID++)
	{
		if ((ITID->second)==1)
		{
			IDA=ITID->first;
		}
		else if ((ITID->second)==2)
		{
			IDB=ITID->first;
		}
	}
	OUT<<"#Chr\tStart\tEnd\tMean_r^2_"<<IDA<<"\tSum_r^2_"<<IDA<<"\tCount_"<<IDA<<"\tMean_r^2_"<<IDB<<"\tSum_r^2_"<<IDB<<"\tCount_"<<IDB<<"\tMeanRRDiff("<<IDA<<"-"<<IDB<<")\n";
	//#####################

	map <string, llong> :: iterator  MaxIT ;
	map <string,map<llong, vector <BaseType>  > > :: iterator key1_it ;

	int  ALL_countV1=0;
	int  ALL_countV2=0;
	double SumMeanRRV1=0.0;
	double SumMeanRRV2=0.0;
	double DiffSum=0.0;
	int ALL_count_Diff=0;

	for (MaxIT=MaxSNPSite.begin(); MaxIT!=MaxSNPSite.end();MaxIT++)
	{
		int MaxBin=int((MaxIT->second)/(paraFA04->bin))+1;
		(paraFA04->count)++;
		StarRsult **All_StatV1 = new StarRsult *[(paraFA04->count)];
		StarRsult **All_StatV2 = new StarRsult *[(paraFA04->count)];
		for (int kk=0 ; kk<(paraFA04->count) ; kk++)
		{
			All_StatV1[kk]=new StarRsult[MaxBin];
			All_StatV2[kk]=new StarRsult[MaxBin];
		}

		key1_it=SNPListV1.find(MaxIT->first);
		PairWiseRRCal( paraFA04,  para_18, key1_it->second ,  All_StatV1);

		key1_it=SNPListV2.find(MaxIT->first);
		PairWiseRRCal( paraFA04,  para_18, key1_it->second ,  All_StatV2);


		MaxBin=MaxBin-(paraFA04->count);
		for (int ii=0 ; ii<MaxBin ; ii++)
		{
			llong Start=ii*(paraFA04->bin);
			llong End=Start+(paraFA04->WindowSize);

			int countV1=0;
			int countV2=0;
			double SumRRV1=0;
			double SumRRV2=0;
			for (int jj=0; jj<(paraFA04->count); jj++)
			{
				for (int uu=jj  ; uu<(paraFA04->count) ; uu++)
				{
					countV1+=(All_StatV1[uu][ii+jj]).Count;
					SumRRV1+=(All_StatV1[uu][ii+jj]).sumRR;
					countV2+=(All_StatV2[uu][ii+jj]).Count;
					SumRRV2+=(All_StatV2[uu][ii+jj]).sumRR;
				}
			}
			double MeanRRV1=-1;
			double MeanRRV2=-1;
			double Diff=-8;			
			if ( (countV1>=(paraFA04->Masked) )   &&  (countV2>=(paraFA04->Masked)) )
			{
				MeanRRV1=SumRRV1/countV1;
				SumMeanRRV1+=MeanRRV1;
				ALL_countV1++;
				MeanRRV2=SumRRV2/countV2;
				SumMeanRRV2+=MeanRRV2;
				ALL_countV2++;
				Diff=MeanRRV1-MeanRRV2 ;
				ALL_count_Diff++;
				DiffSum+=Diff;
			}
			else if  ( (countV1>=(paraFA04->Masked))  &&  (countV2<(paraFA04->Masked))  )
			{
				MeanRRV1=SumRRV1/countV1;
				SumMeanRRV1+=MeanRRV1;
				ALL_countV1++;
			}
			else if ( (countV1<(paraFA04->Masked))  &&  (countV2>=(paraFA04->Masked))  )
			{
				MeanRRV2=SumRRV2/countV2;
				SumMeanRRV2+=MeanRRV2;
				ALL_countV2++;
			}

			OUT<<MaxIT->first<<"\t"<<Start<<"\t"<<End<<"\t"<<setprecision(4)<<setiosflags(ios::left)<<setiosflags(ios::fixed)<<MeanRRV1<<"\t"<<SumRRV1<<"\t"<<countV1<<"\t"<<MeanRRV2<<"\t"<<SumRRV2<<"\t"<<countV2<<"\t"<<Diff<<"\n";
		}
		for (int kk=0 ; kk<(paraFA04->count) ; kk++)
		{
			delete [] All_StatV1[kk];
			delete [] All_StatV2[kk];
		}
		delete [] All_StatV1 ;
		delete [] All_StatV2 ;
	}
	OUT.close();



	double meaWolRRV1=SumMeanRRV1/ALL_countV1;
	double meaWolRRV2=SumMeanRRV2/ALL_countV2;
	double meaWolDiff=DiffSum/ALL_count_Diff;

	igzstream INAA (Stat.c_str(),ifstream::in);
	if (INAA.fail())
	{
		cerr << "open Sub Group IN File error: "<<Stat<<endl;
		return  0;
	}


	string StatSS=(paraFA04->InStr2)+".winRR.gz";
	ogzstream OUTSS ((StatSS).c_str());
	if((!OUTSS.good()))
	{
		cerr << "open OUT File error: "<<Stat<<endl;
		return 0;
	}

	string StatSELE=(paraFA04->InStr2)+".diffRegion.gz";
	ogzstream OUTSELE ((StatSELE).c_str());
	if((!OUTSELE.good()))
	{
		cerr << "open OUT File error: "<<StatSELE<<endl;
		return 0;
	}



	string chr;
	string start;
	string end;
	getline(INAA,chr);

	long double diffV1=0;
	long double diffV2=0;
	long double diffDiff=0;
	long double SDV1=0;
	long double SDV2=0;
	long double SDDiff=0;
	double AARRV1;
	double AARRV2;

	string CountV1, SumESV1 ;
	string CountV2, SumESV2 ;
	double DiffHere ;
	//OUT<<"#Chr\tStart\tEnd\tMean_r^2_"<<IDA<<"\tSum_r^2"<<IDA<<"\tCount"<<IDA<<"\tMean_r^2_"<<IDB<<"\tSum_r^2"<<IDB<<"\tCount"<<IDB<<"\n";
	while(!INAA.eof())
	{
		string  line ;
		getline(INAA,line);
		if (line.length()<=0)  { continue ; }
		istringstream isone (line,istringstream::in);
		isone>>chr>>start>>end>>AARRV1>>SumESV1>>CountV1>>AARRV2>>SumESV2>>CountV2>>DiffHere;
		if (AARRV1>=0 )  
		{
			diffV1=AARRV1-meaWolRRV1;
			SDV1+=(diffV1*diffV1);
		}
		if (AARRV2>=0) 
		{
			diffV2=AARRV2-meaWolRRV2;
			SDV2+=(diffV2*diffV2);
		}
		if ( DiffHere>-6)
		{
			diffDiff=DiffHere-meaWolDiff;
			SDDiff+=(diffDiff*diffDiff);
		}

	}
	INAA.close();
	SDV1=sqrt(SDV1/ALL_countV1);
	SDV2=sqrt(SDV2/ALL_countV2);
	SDDiff=sqrt(SDDiff/ALL_count_Diff);


	igzstream INBB (Stat.c_str(),ifstream::in);
	if (INBB.fail())
	{
		cerr << "open Sub Group IN File error: "<<Stat<<endl;
		return  0;
	}

	getline(INBB,chr);
	double PValue, ZScore;
	OUTSS<<chr<<"\tZScore\tPvalue\t\n";
	OUTSELE<<chr<<"\tZScore\tPvalue\t\n";
	OUTSS<<"##Group["<<IDA<<"], MeanRR:"<<meaWolRRV1<<"\tSD:"<<SDV1<<"\tEffective windows Count:"<<ALL_countV1<<"\n";
	OUTSS<<"##Group["<<IDB<<"], MeanRR:"<<meaWolRRV2<<"\tSD:"<<SDV2<<"\tEffective windows Count:"<<ALL_countV2<<"\n";
	OUTSS<<"##Diff MeanRR["<<IDA<<"-"<<IDB<<"], Mean:"<<meaWolDiff<<"\tSD:"<<SDDiff<<"\tEffective windows Count:"<<ALL_count_Diff<<endl;

	while(!INBB.eof())
	{
		string  line ;
		getline(INBB,line);
		if (line.length()<=0) { continue ;}
		istringstream isone (line,istringstream::in);
		isone>>chr>>start>>end>>AARRV1>>SumESV1>>CountV1>>AARRV2>>SumESV2>>CountV2>>DiffHere;
		if  ( DiffHere> -6 )
		{ 
			diffV1=DiffHere-meaWolDiff;
			PValue=(1-pnorm5(DiffHere, meaWolDiff , SDDiff , 1 , 0 ));
			ZScore=diffV1/SDDiff;
			OUTSS<<line<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::left)<<setprecision(2)<<ZScore<<"\t"<<setprecision(6)<<PValue<<endl;
			if (PValue<(paraFA04->Pvalue))
			{
				OUTSELE<<line<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::left)<<setprecision(2)<<ZScore<<"\t"<<setprecision(6)<<PValue<<endl;
			}
		}
		else
		{
			if  (AARRV1 <0  &&  AARRV2 < 0 )
			{
				OUTSS<<chr<<"\t"<<start<<"\t"<<end<<"\tNA\t"<<SumESV1<<"\t"<<CountV1<<"\tNA\t"<<SumESV2<<"\t"<<CountV2<<"\tNA\n";
			}
			else if ( AARRV1 <0   &&  AARRV2>=0)
			{
				OUTSS<<chr<<"\t"<<start<<"\t"<<end<<"\tNA\t"<<SumESV1<<"\t"<<CountV1<<"\t"<<AARRV2<<"\t"<<SumESV2<<"\t"<<CountV2<<"\tNA\n";				
			}
			else
			{
				OUTSS<<chr<<"\t"<<start<<"\t"<<end<<"\t"<<AARRV1<<"\t"<<SumESV1<<"\t"<<CountV1<<"\tNA\t"<<SumESV2<<"\t"<<CountV2<<"\tNA\n";
			}
		}
	}
	INBB.close();
	OUTSS.close();
	OUTSELE.close();

	string cc="rm  "+Stat ;
	std::system(cc.c_str()) ;

	DrawTowDepthSliding(StatSS, paraFA04, para_18 , IDA , IDB);
	return 1;

}







//////////////////////////////////////////////////////////////////  one Group  Deal ////////////////////////////////////////////////////
//
//
//
//
//
//
//
//
//
//////////////////////////////////////////////////////////////////  one Group  Deal ///////////////////////////////////////////////////






void DrawDepthSliding( string File , In3str1v * paraFA04, Para_18 * para_18 )
{

	//  OUTSS<<"Chr\tStart\tEnd\tMeanRR\tSum_r^2\tCount\tZScore\tPvalue\n";

	string OutPlotTmp=(paraFA04->InStr2)+".tmp";
	ofstream  OUTPLOT (OutPlotTmp.c_str());

	igzstream LIST (File.c_str(),ifstream::in);

	string chr="";
	getline(LIST,chr);
	getline(LIST,chr);
	string Start,End,MeanRR;
	map <string , int >  VecChr;
	map <string , map <int , string > >  Value;
	map <string , map <int , string > > :: iterator iTValue;

	map  <string , int >  :: iterator itVecChr ;
	while(!LIST.eof())
	{
		string  line ;
		getline(LIST,line);
		if (line.length()<=0) { continue ;}
		istringstream isone (line,istringstream::in);
		isone>>chr>>Start>>End>>MeanRR ;
		itVecChr=VecChr.find(chr);
		if (itVecChr==VecChr.end())
		{
			VecChr.insert(map <string,int>::value_type(chr,1));
			map <int , string >  DataEE;
			DataEE[0]=MeanRR;
			Value.insert(map <string , map <int , string > > ::value_type(chr,DataEE));
		}
		else
		{
			iTValue=Value.find(chr);
			(iTValue->second).insert(map <int , string >::value_type(itVecChr->second,MeanRR));
			itVecChr->second++;
		}
	}
	LIST.close();

	int MaxXX=0;
	OUTPLOT<<"Site";
	for(itVecChr=VecChr.begin();itVecChr!=VecChr.end();itVecChr++)
	{
		if (MaxXX<(itVecChr->second)) {MaxXX=itVecChr->second;}
		OUTPLOT<<"\t"<<itVecChr->first;
	}
	OUTPLOT<<endl;


	map <int , string > ::  iterator  itEEEE;
	for (int ii=0 ; ii< MaxXX ; ii++)
	{
		OUTPLOT<<ii;
		for(itVecChr=VecChr.begin();itVecChr!=VecChr.end();itVecChr++)
		{
			iTValue=Value.find(itVecChr->first);
			itEEEE=(iTValue->second).find(ii);

			if (itEEEE==(iTValue->second).end()) 
			{
				OUTPLOT<<"\tNA";
			}
			else
			{
				OUTPLOT<<"\t"<<itEEEE->second ;
			}
		}
		OUTPLOT<<endl;

	}


	OUTPLOT.close();

	string OutPlotr=(paraFA04->InStr2)+".r";
	ofstream  OUTR (OutPlotr.c_str());

	OUTR<<""
		"\n"
		"library(ggplot2);\nlibrary(reshape2);\n"
		"read.table(\""<<OutPlotTmp<<"\",header=T)->data;\n"
		"data2<-melt(data,id=1);\n"
		"p <- ggplot(data = data2)+geom_tile(aes(x = Site,y = variable, fill=value))+scale_fill_gradient(limit=c(0,1),low ='white',high ='red',na.value=\"grey50\")\n"
		"pdf(\""<<(paraFA04->InStr2)<<".pdf\",h=6,w=25);\n"
		"p\n"
		"dev.off()\n"
		"png(\""<<(paraFA04->InStr2)<<".png\",h=6,w=25);\n"
		"p\n"
		"dev.off()\n"<<endl;
	OUTR.close();

	char   buf[2048]={'\0'};
	string cc="which  Rscript  2> /dev/null ";
	memset( buf, '\0', sizeof(buf) );
	FILE   *stream ;
	stream=popen(cc.c_str(),"r") ;
	fread( buf, sizeof(char), sizeof(buf), stream);
	string binPath=buf;
	binPath=binPath.substr(0,binPath.length()-1);
	if (binPath == "" )
	{
		cout <<"\twarning: can't find the [Rscript] in your $PATH ; no png Figure Out"<<endl;
		cout <<"\t\tRscript "<<OutPlotr<<endl;
	}
	else
	{
		if (paraFA04->TF)
		{
			//cc=binPath+"\t"+OutPlotr ;
			cc=binPath+"\t"+OutPlotr+"  ; rm  -rf  "+ OutPlotr ;
			std::system(cc.c_str()) ;
		}
		else
		{
			cc=binPath+"\t"+OutPlotr ;
			std::system(cc.c_str()) ;
			cout <<"\t\tRePlot : Rscript  "<<OutPlotr<<endl;
		}
	}
}






int SilidingRRCal( In3str1v * paraFA04, Para_18 * para_18 , map <string,map <llong, vector <BaseType> > >  & SNPList ,  map <string,llong>  & MaxSNPSite)
{
	string Stat=(paraFA04->InStr2)+".winRR.temp.gz";
	ogzstream OUT ((Stat).c_str());
	if((!OUT.good()))
	{
		cerr << "open OUT File error: "<<Stat<<endl;
		return 0;
	}
	OUT<<"#Chr\tStart\tEnd\tMean_r^2\tSum_r^2\tCount\n";
	//#####################

	map <string, llong> :: iterator  MaxIT ;
	map <string,map<llong, vector <BaseType>  > > :: iterator key1_it ;

	int  ALL_count=0;
	double SumMeanRR=0.0;

	for (MaxIT=MaxSNPSite.begin(); MaxIT!=MaxSNPSite.end();MaxIT++)
	{
		int MaxBin=int((MaxIT->second)/(paraFA04->bin))+1;
		(paraFA04->count)++;
		StarRsult **All_Stat = new StarRsult *[(paraFA04->count)];
		for (int kk=0 ; kk<(paraFA04->count) ; kk++)
		{
			All_Stat[kk]=new StarRsult[MaxBin];
		}

		key1_it=SNPList.find(MaxIT->first);
		PairWiseRRCal( paraFA04,  para_18, key1_it->second ,  All_Stat);

		MaxBin=MaxBin-(paraFA04->count);
		for (int ii=0 ; ii<MaxBin ; ii++)
		{
			int count=0;
			double SumRR=0;
			for (int jj=0; jj<(paraFA04->count); jj++)
			{
				for (int uu=jj  ; uu<(paraFA04->count) ; uu++)
				{
					count+=(All_Stat[uu][ii+jj]).Count;
					SumRR+=(All_Stat[uu][ii+jj]).sumRR;
				}
			}
			if (count<(paraFA04->Masked))  {continue;}
			double MeanRR=SumRR/count;
			llong Start=ii*(paraFA04->bin);
			llong End=Start+(paraFA04->WindowSize);
			SumMeanRR+=MeanRR;
			ALL_count++;
			OUT<<MaxIT->first<<"\t"<<Start<<"\t"<<End<<"\t"<<setprecision(4)<<setiosflags(ios::left)<<setiosflags(ios::fixed)<<MeanRR<<"\t"<<SumRR<<"\t"<<count<<"\n";
		}

		for (int kk=0 ; kk<(paraFA04->count) ; kk++)
		{
			delete [] All_Stat[kk];
		}

		delete [] All_Stat ;
	}
	OUT.close();



	double meaWolRR=SumMeanRR/ALL_count;


	igzstream INAA (Stat.c_str(),ifstream::in);
	if (INAA.fail())
	{
		cerr << "open Sub Group IN File error: "<<Stat<<endl;
		return  0;
	}


	string StatSS=(paraFA04->InStr2)+".winRR.gz";
	ogzstream OUTSS ((StatSS).c_str());
	if((!OUTSS.good()))
	{
		cerr << "open OUT File error: "<<Stat<<endl;
		return 0;
	}

	string StatSELE=(paraFA04->InStr2)+".selection.gz";
	ogzstream OUTSELE ((StatSELE).c_str());
	if((!OUTSELE.good()))
	{
		cerr << "open OUT File error: "<<StatSELE<<endl;
		return 0;
	}



	double AARR;
	string chr;
	string start;
	string end;
	getline(INAA,chr);

	long double diff=0;
	long double SD=0;

	while(!INAA.eof())
	{
		string  line ;
		getline(INAA,line);
		if (line.length()<=0)  { continue ; }
		istringstream isone (line,istringstream::in);
		isone>>chr>>start>>end>>AARR;
		diff=AARR-meaWolRR;
		SD+=(diff*diff);
	}
	INAA.close();
	SD=sqrt(SD/ALL_count);


	igzstream INBB (Stat.c_str(),ifstream::in);
	if (INBB.fail())
	{
		cerr << "open Sub Group IN File error: "<<Stat<<endl;
		return  0;
	}

	getline(INBB,chr);
	double PValue, ZScore;
	OUTSS<<"#Chr\tStart\tEnd\tMeanRR\tSum_r^2\tCount\tZScore\tPvalue\n";
	OUTSELE<<"#Chr\tStart\tEnd\tMean_r^2\tSum_r^2\tCount\tZScore\tPvalue\n";
	OUTSS<<"##MeanRR:"<<meaWolRR<<"\tSD:"<<SD<<"\tEffective windows Count:"<<ALL_count<<endl;

	while(!INBB.eof())
	{
		string  line ;
		getline(INBB,line);
		if (line.length()<=0) { continue ;}
		istringstream isone (line,istringstream::in);
		isone>>chr>>start>>end>>AARR;
		diff=AARR-meaWolRR;
		PValue=(1-pnorm5(AARR, meaWolRR , SD  , 1 , 0 ));
		ZScore=diff/SD;
		if  (ZScore<0){ZScore=0-ZScore;}
		OUTSS<<line<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::left)<<setprecision(2)<<ZScore<<"\t"<<setprecision(6)<<PValue<<endl;
		if  ( diff>0  &&  PValue<(paraFA04->Pvalue))
		{
			OUTSELE<<line<<"\t"<<setiosflags(ios::fixed)<<setiosflags(ios::left)<<setprecision(2)<<ZScore<<"\t"<<setprecision(6)<<PValue<<endl;
		}
	}
	INBB.close();
	OUTSS.close();
	OUTSELE.close();

	string cc="rm  "+Stat ;
	std::system(cc.c_str()) ;

	DrawDepthSliding(StatSS, paraFA04, para_18 );
	return 1;

}





#endif // processV3_H_ //
///////// swimming in the sky and flying in the sea ////////////




