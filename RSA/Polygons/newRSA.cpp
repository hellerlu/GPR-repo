#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <ctime>
#include <cstddef>
#include <memory>
#include <cstring>
#define NDEBUG
#include <cassert>


#include "RandomGenerator.h"


const int MaxDimension=4;
const double SystemSize=1;



size_t TimeLimit;
time_t start;

const double Pi = 3.1415926535897932384626433832795028841971693993751;
double UnitSphereVolume(int dimension)//return the volume of a hypersphere of radius 1
{
	if (dimension == 1) return 2;
	else if (dimension == 2) return Pi;
	else if (dimension == 3) return 4 * Pi / 3;
	else if (dimension == 4) return Pi*Pi / 2;
	else if (dimension == 5) return 8 * Pi*Pi / 15;
	else if (dimension == 6) return Pi*Pi*Pi / 6;
	else if (dimension == 7) return 16 * Pi*Pi*Pi / 105;
	else if (dimension == 8) return Pi*Pi*Pi*Pi / 24;
	else
		throw "error in SphereVolume: unsupported dimension";
}

class progress_display
{
public:
	explicit progress_display(unsigned long expected_count,
		std::ostream & os = std::cout,
		const std::string & s1 = "\n", //leading strings
		const std::string & s2 = "",
		const std::string & s3 = "")
		// os is hint; implementation may ignore, particularly in embedded systems
		: m_os(os), m_s1(s1), m_s2(s2), m_s3(s3) {
		restart(expected_count);
	}

	void restart(unsigned long expected_count)
		//  Effects: display appropriate scale
		//  Postconditions: count()==0, expected_count()==expected_count
	{
		_count = _next_tic_count = _tic = 0;
		_expected_count = expected_count;

		m_os << m_s1 << "0%   10   20   30   40   50   60   70   80   90   100%\n"
			<< m_s2 << "|----|----|----|----|----|----|----|----|----|----|"
			<< std::endl  // endl implies flush, which ensures display
			<< m_s3;
		if (!_expected_count) _expected_count = 1;  // prevent divide by zero
	} // restart

	unsigned long  operator+=(unsigned long increment)
		//  Effects: Display appropriate progress tic if needed.
		//  Postconditions: count()== original count() + increment
		//  Returns: count().
	{
		if ((_count += increment) >= _next_tic_count) { display_tic(); }
		return _count;
	}

	unsigned long  operator++()           { return operator+=(1); }
	unsigned long  operator++(int)           { return (operator+=(1) - 1); }
	unsigned long  count() const          { return _count; }
	unsigned long  expected_count() const { return _expected_count; }

private:
	std::ostream &     m_os;  // may not be present in all imps
	const std::string  m_s1;  // string is more general, safer than 
	const std::string  m_s2;  //  const char *, and efficiency or size are
	const std::string  m_s3;  //  not issues

	unsigned long _count, _expected_count, _next_tic_count;
	unsigned int  _tic;
	void display_tic()
	{
		// use of floating point ensures that both large and small counts
		// work correctly.  static_cast<>() is also used several places
		// to suppress spurious compiler warnings. 
		unsigned int tics_needed =
			static_cast<unsigned int>(
			(static_cast<double>(_count) / _expected_count)*50.0);
		do { m_os << '*' << std::flush; } while (++_tic < tics_needed);
		_next_tic_count =
			static_cast<unsigned long>((_tic / 50.0)*_expected_count);
		if (_count == _expected_count) {
			if (_tic < 51) m_os << '*';
			m_os << std::endl;
		}
	} // display_tic
};


class CellList;
class VoxelList;
class Shape
{
public:
	//dimension of the shape
	virtual int GetDimension() = 0;

	//number of parameters to uniquely define a particle of the shape
	//for fixed shapes, this is usually equal to dimension+(rotational degrees of freedom)
	virtual int GetVoxelDimension() = 0;

	//volume of a particle
	virtual double GetVolume() = 0;

	//if two particles are farther than this distance, they cannot overlap
	virtual double GetInteractionRadius() = 0;

	//test if two particles overlap
	virtual bool Overlap(const double * adjParticle, const double * Coordinate) = 0;

	//test if a voxel at Coordinate is completely covered by adjParticle
	virtual bool VoxelOverlap(const double * adjParticle, const double * Coordinate, double HalfSize, double orientHalfSize) = 0;

	//some shapes require additional processing after voxels are splitted, but the default is an empty function
	virtual void ProcessSplittedVoxels(const CellList & Configuration, VoxelList & voxels, std::ostream & output, RandomGenerator & gen)
	{
	}

	virtual ~Shape(){}
};

struct CellList_TempStruct
{
	signed char c[MaxDimension];
	float t;
};
bool operator < (const CellList_TempStruct & Left, const CellList_TempStruct & Right)
{
	return Left.t<Right.t;
}
class CellList
{
private:
	int dimension;
	signed int cellrank;
	unsigned int nbradjacentcell;
	double cellsize;//side length of a cell;
	unsigned int * head;//head pointer of the static linklist
	unsigned int * next;//static linklist
	signed char * adjacentcelllist;
public:
	unsigned int nbrparticle, maxnbrparticle;
	double * coordinate;//static linklist
	unsigned int getindex2(const double * coord) const//find index of the cell corresponding to coord
	{
		signed int index[MaxDimension];
		this->getindex(coord, index);
		return this->convertindex(index);
	}
	unsigned int convertindex(const signed int * index) const//convert dimensional index to the 1-dimension cell index
	{
		unsigned int result=index[0];
		for(int i=1; i<this->dimension; i++)
		{
			result*=this->cellrank;
			result+=index[i];
		}
		return result;
	}
	void getindex(const double * coord, signed int * index) const
	{
		for(int i=0; i<this->dimension; i++)
			index[i]=std::floor(coord[i]/this->cellsize);
	}
	void getadjacentcelllist(void)
	{
		int max=(int)(std::pow((double)(3),this->dimension)+0.1);
		this->nbradjacentcell=max;
		this->adjacentcelllist= new signed char [max*MaxDimension];
		signed char increment=1;
		signed char finish=1;
	    signed char coord[MaxDimension];
		float temp;
	    int i;
		std::vector<CellList_TempStruct> result;
		result.resize(max);
	    std::vector<CellList_TempStruct>::iterator iter=result.begin();
	    //initialize
		for(i=0; i<this->dimension; i++)
			coord[i]=-1;
	    //main loop body
	start:
		temp=0;
	    for(i=0; i<this->dimension; i++)
		{
	        iter->c[i]=coord[i];
			temp+=coord[i]*coord[i];
		}
		iter->t=temp;	
	    iter++;
	    //loop end
	    coord[0]+=increment;
	    for(i=0; i<this->dimension-1; i++)
	        if(coord[i]>finish)
	        {
	            coord[i]=-1;
	            coord[i+1]+=increment;
	        }
	    if(coord[this->dimension-1]<=finish) goto start;
		// end of the loop

		std::sort(result.begin(),result.end());
		//copy the data from vector result to array adjacentcelllist
		for(i=0; i<max; i++)
		{
			signed char * now=this->adjacentcelllist+i*MaxDimension;
			int j;
			for(j=0; j<this->dimension; j++)
				now[j]=result[i].c[j];
			for(; j<MaxDimension; j++)
				now[j]=0;
		}

	    return;
	}



	CellList()
	{
		this->dimension=0;
		this->nbrparticle=0;
		this->maxnbrparticle=0;
		this->head=nullptr;
		this->next=nullptr;
		this->coordinate=nullptr;
		this->adjacentcelllist=nullptr;
	}
	CellList(int Dimension, unsigned int MaxNbrSphere, double InteractionRadius)
	{
		this->dimension=Dimension;
		if(this->dimension>MaxDimension) throw "Error initializing celllist: dimension greater than max dimension";
		this->nbrparticle=0;
		this->maxnbrparticle=MaxNbrSphere;
		this->cellrank=std::floor(SystemSize/InteractionRadius);
		this->cellsize=SystemSize/this->cellrank;
		size_t temp=1;
		for(int i=0; i<this->dimension; i++)
			temp*=this->cellrank;
		this->head=new unsigned int[temp];
		for(size_t i=0; i<temp; i++)
			this->head[i]=0;//0 is the end of static linklist
		this->next=new unsigned int[this->maxnbrparticle];
		this->coordinate=new double[MaxDimension*this->maxnbrparticle];
		this->getadjacentcelllist();
	}	
	CellList(const CellList & src)
	{
		std::memcpy(this, &src, sizeof(CellList));

		size_t temp=1;
		for(int i=0; i<this->dimension; i++)
			temp*=this->cellrank;
		this->head=new unsigned int[temp];
		std::memcpy(this->head, src.head, sizeof(unsigned int)*temp);
		this->next=new unsigned int[this->maxnbrparticle];
		std::memcpy(this->next, src.next, sizeof(unsigned int)*src.maxnbrparticle);
		this->coordinate=new double[MaxDimension*this->maxnbrparticle];
		std::memcpy(this->coordinate, src.coordinate, sizeof(double)*::MaxDimension*src.maxnbrparticle);
		this->adjacentcelllist=nullptr;
		this->getadjacentcelllist();
	}
	~CellList()
	{
		if(this->head!=nullptr)delete[] this->head;
		if(this->next!=nullptr)delete[] this->next;
		if(this->coordinate!=nullptr)delete[] this->coordinate;
		if(this->adjacentcelllist!=nullptr)delete[] this->adjacentcelllist;
	}
	void Insert(const double * Coordinate)
	{
		this->nbrparticle++;//ignore the 0th particle since 0 is used as a flag
		if(this->nbrparticle>=this->maxnbrparticle) throw "Error in celllist: out of capacity";

		double * data=this->coordinate+MaxDimension*this->nbrparticle;

		//copy coordinate
		for(int i=0; i<MaxDimension; i++)//copy coordinate
			data[i]=Coordinate[i];

		//insert linklist
		unsigned int index=this->getindex2(Coordinate);
		this->next[this->nbrparticle]=this->head[index];
		this->head[index]=this->nbrparticle;
	}
	void FindNeighbour(const double * Coordinate) const //iterate through neighbours in adjacent cells, do nothing
	{
		signed int index[MaxDimension], index2[MaxDimension];
		signed char out[MaxDimension];//if the coordinate is out of system and in its periodic image
		double adjparticle[MaxDimension];//the coordinate of this adjacent particle
		if(this->cellrank>3)
		{
			this->getindex(Coordinate, index);
			for(unsigned int i=0; i<this->nbradjacentcell; i++)
			{
				signed char * now=this->adjacentcelllist+MaxDimension*i;
				for(int j=0; j<this->dimension; j++)
				{
					//generate index2:the index of current adjacentcell
					index2[j]=index[j]+now[j];
					if(index2[j]<0)
					{
						index2[j]+=this->cellrank;
						out[j]=-1;
					}
					else if(index2[j]>=this->cellrank)
					{
						index2[j]-=this->cellrank;
						out[j]=1;
					}
					else
						out[j]=0;
				}

				unsigned int cparticle;//current particle's index in the static linklist
				cparticle=head[this->convertindex(index2)];
				while(cparticle!=0)//0 is the flag of end of a linklist
				{
					//generate coordinate of the new particle
					for(int j=0; j<this->dimension; j++)
						adjparticle[j]=this->coordinate[cparticle*MaxDimension+j]+out[j]*SystemSize;
						
					//add whatever you want to do here



					cparticle=this->next[cparticle];
				}
			}
			return ;
		}
		else
		{
			//cell list is not helping
			for(unsigned int i=1; i<this->nbrparticle+1; i++)
			{
				for(int j=0; j<this->dimension; j++)
				{
					adjparticle[j]=this->coordinate[i*MaxDimension+j];
					if((adjparticle[j] - Coordinate[j]) < (-0.5)*::SystemSize)
						adjparticle[j]+=::SystemSize;
					else if((adjparticle[j] - Coordinate[j]) > (0.5)*::SystemSize)
						adjparticle[j]-=::SystemSize;
				}

				//add whatever you want to do here


			}
			return ;
		}
	}
	bool CheckOverlap(const double * Coordinate, Shape * pShape) const
	{
		signed int index[MaxDimension], index2[MaxDimension];
		signed char out[MaxDimension] = { 0 };//if the coordinate is out of system and in its periodic image
		double adjparticle[MaxDimension];//the coordinate of this adjacent particle
		if(this->cellrank>3)
		{
			this->getindex(Coordinate, index);
			for(unsigned int i=0; i<this->nbradjacentcell; i++)
			{
				signed char * now=this->adjacentcelllist+MaxDimension*i;
				for(int j=0; j<this->dimension; j++)
				{
					//generate index2:the index of current adjacentcell
					index2[j]=index[j]+now[j];
					if(index2[j]<0)
					{
						index2[j]+=this->cellrank;
						out[j]=-1;
					}
					else if(index2[j]>=this->cellrank)
					{
						index2[j]-=this->cellrank;
						out[j]=1;
					}
					else
						out[j]=0;
				}

				unsigned int cparticle;//current particle's index in the static linklist
				cparticle=head[this->convertindex(index2)];
				while(cparticle!=0)//0 is the flag of end of a linklist
				{
					//generate coordinate of the new particle
					for(int j=0; j<MaxDimension; j++)
						adjparticle[j]=this->coordinate[cparticle*MaxDimension+j]+out[j]*SystemSize;
						
					//add whatever you want to do here
					if (pShape->Overlap(adjparticle, Coordinate))
						return true;

					cparticle=this->next[cparticle];
				}
			}
			return false;
		}
		else
		{
			//cell list is not helping
			for(unsigned int i=1; i<this->nbrparticle+1; i++)
			{
				for(int j=0; j<MaxDimension; j++)
				{
					adjparticle[j]=this->coordinate[i*MaxDimension+j];
					if((adjparticle[j] - Coordinate[j]) < (-0.5)*::SystemSize)
						adjparticle[j]+=::SystemSize;
					else if((adjparticle[j] - Coordinate[j]) > (0.5)*::SystemSize)
						adjparticle[j]-=::SystemSize;
				}

				//add whatever you want to do here
				if (pShape->Overlap(adjparticle, Coordinate))
					return true;
			}
			return false;
		}
	}
	bool CheckVoxelOverlap(const double * Coordinate, Shape * pShape, double HalfSize, double orientHalfSize) const//Check for voxel overlap, Coordinate is the coordinate of the voxel CENTER
	{
		signed int index[MaxDimension], index2[MaxDimension];
		signed char out[MaxDimension] = { 0 };//if the coordinate is out of system and in its periodic image
		double adjparticle[MaxDimension];//the coordinate of this adjacent particle
		if(this->cellrank>3)
		{
			this->getindex(Coordinate, index);
			for(unsigned int i=0; i<this->nbradjacentcell; i++)
			{
				signed char * now=this->adjacentcelllist+MaxDimension*i;
				for(int j=0; j<this->dimension; j++)
				{
					//generate index2:the index of current adjacentcell
					index2[j]=index[j]+now[j];
					if(index2[j]<0)
					{
						index2[j]+=this->cellrank;
						out[j]=-1;
					}
					else if(index2[j]>=this->cellrank)
					{
						index2[j]-=this->cellrank;
						out[j]=1;
					}
					else
						out[j]=0;
				}

				unsigned int cparticle;//current particle's index in the static linklist
				cparticle=head[this->convertindex(index2)];
				while(cparticle!=0)//0 is the flag of end of a linklist
				{
					//generate coordinate of the new particle
					for(int j=0; j<MaxDimension; j++)
						adjparticle[j]=this->coordinate[cparticle*MaxDimension+j]+out[j]*SystemSize;
						
					//add whatever you want to do here
					if (pShape->VoxelOverlap(adjparticle, Coordinate, HalfSize, orientHalfSize))
						return true;

					cparticle=this->next[cparticle];
				}
			}
			return false;
		}
		else
		{
			//cell list is not helping
			for(unsigned int i=1; i<this->nbrparticle+1; i++)
			{
				for(int j=0; j<MaxDimension; j++)
				{
					adjparticle[j]=this->coordinate[i*MaxDimension+j];
					if((adjparticle[j] - Coordinate[j]) < (-0.5)*::SystemSize)
						adjparticle[j]+=::SystemSize;
					else if((adjparticle[j] - Coordinate[j]) > (0.5)*::SystemSize)
						adjparticle[j]-=::SystemSize;
				}

				//add whatever you want to do here
				if (pShape->VoxelOverlap(adjparticle, Coordinate, HalfSize, orientHalfSize))
					return true;
			}
			return false;
		}
	}


	template<typename OutType> void OutputCoordinate(OutType & out) const
	{
		for(unsigned int i=1; i<this->nbrparticle+1; i++)
		{
			double * now = this->coordinate+i*MaxDimension;
			for(int j=0; j<this->dimension; j++)
				out<<now[j]<<'\t';
			out<<'\n';
		}
	}


	void Save(std::ostream & writeto) const
	{
		writeto.write((char*)(&this->nbrparticle), sizeof(this->nbrparticle));
		for(unsigned int i=1; i<this->nbrparticle+1; i++)
		{
			double * now = this->coordinate+i*MaxDimension;
			writeto.write((char*)(now), sizeof(double)*MaxDimension);
		}
	}
	void Load(std::istream & source)
	{
		unsigned int nbr;
		source.read((char*)(&nbr), sizeof(this->nbrparticle));
		double temp[MaxDimension];
		for(unsigned int i=0; i<nbr; i++)
		{
			source.read((char*)(temp), sizeof(double)*MaxDimension);
			this->Insert(temp);
		}
	}


};

template <typename T> void SaveVectorData(const std::vector<T> & data, std::ostream & writeto)
{
	size_t a=data.size();
	writeto.write((char*)(&a), sizeof(a));
	writeto.write((char*)(&data[0]), a*sizeof(T));
}
template <typename T> void LoadVectorData(std::vector<T> & data, std::istream & readfrom)
{
	size_t a;
	readfrom.read((char*)(&a), sizeof(a));
	data.resize(a);
	readfrom.read((char*)(&data[0]), a*sizeof(T));
}

#include <omp.h>
typedef unsigned short voxeltype;
class VoxelList
{
private:
	int ParticleDimension, VoxelDimension, level;
	double halfsize, originhalfsize;
	double orienthalfsize, orientoriginhalfsize;
	std::vector<voxeltype> * originindex;
	std::vector<unsigned short> * subindex;
	unsigned short nbrsplitvoxel;
	signed char * splitvoxellist;
	void getsplitvoxellist(void)
	{
		int max=(int)(std::pow((double)(2),this->VoxelDimension)+0.1);
		this->nbrsplitvoxel=max;
		this->splitvoxellist= new signed char [max*MaxDimension];
		signed char increment=2;
		signed char finish=1;
	    signed char coord[MaxDimension];
	    int i;
		std::vector<CellList_TempStruct> result;
		result.resize(max);
	    std::vector<CellList_TempStruct>::iterator iter=result.begin();
	    //initialize
		for (i = 0; i<this->VoxelDimension; i++)
			coord[i]=-1;
	    //main loop body
	start:
		for (i = 0; i<this->VoxelDimension; i++)
		{
	        iter->c[i]=coord[i];
		}
	    iter++;
	    //loop end
	    coord[0]+=increment;
		for (i = 0; i<this->VoxelDimension - 1; i++)
	        if(coord[i]>finish)
	        {
	            coord[i]=-1;
	            coord[i+1]+=increment;
	        }
		if (coord[this->VoxelDimension - 1] <= finish) goto start;
		// end of the loop

		//copy the data from vector result to array adjacentcelllist
		for(i=0; i<max; i++)
		{
			signed char * now=this->splitvoxellist+i*MaxDimension;
			int j;
			for (j = 0; j<this->VoxelDimension; j++)
				now[j]=result[i].c[j];
			for(; j<MaxDimension; j++)
				now[j]=0;
		}

	    //function end
	    return;
	}

public:
	VoxelList()
	{
		this->ParticleDimension = 0;
		this->VoxelDimension = 0;
		this->level=0;
		this->halfsize=0;
		this->originhalfsize=0;
		this->nbrsplitvoxel=0;
		this->originindex=nullptr;
		this->subindex=nullptr;
		this->splitvoxellist=nullptr;
	}
	~VoxelList()
	{
		if(this->originindex!=nullptr) delete this->originindex;
		if(this->subindex!=nullptr) delete this->subindex;
		if(this->splitvoxellist!=nullptr) delete [] this->splitvoxellist;
	}
	VoxelList(int ParticleDimension, int VoxelDimension, voxeltype VoxelRank, voxeltype orientVoxelRank, size_t ExpectedSize=10000)
		: ParticleDimension(ParticleDimension), VoxelDimension(VoxelDimension)
	{
		this->level=0;
		this->halfsize=::SystemSize/VoxelRank/2;
		this->originhalfsize=this->halfsize;
		this->orienthalfsize = 1.0 / orientVoxelRank / 2;
		this->orientoriginhalfsize = this->orienthalfsize;
		this->getsplitvoxellist();
		this->originindex=new std::vector<voxeltype>;
		this->originindex->reserve(ExpectedSize*MaxDimension);
		this->subindex=new std::vector<unsigned short>;
	}
	void GetOriginalListSerial(const CellList & Configuration, Shape * pShape)
	{
		short i;
		double increment[MaxDimension];
		for(i=0; i<this->ParticleDimension; i++)
			increment[i]= 2.0*this->halfsize;
		for(; i<this->VoxelDimension; i++)
			increment[i]= 2.0*this->orienthalfsize;
		double finish = ::SystemSize;
		double initial[MaxDimension];
		for (i = 0; i < this->ParticleDimension; i++)
			initial[i] = this->halfsize;
		for (; i < this->VoxelDimension; i++)
			initial[i] = this->orienthalfsize;
		double coord[MaxDimension];
		//initialize
		for (i = 0; i<this->VoxelDimension; i++)
			coord[i] = initial[i];
		//main loop body
	start:
		if (Configuration.CheckVoxelOverlap(coord, pShape, this->halfsize, this->orienthalfsize) == false)
		{
			for (i = 0; i<this->VoxelDimension; i++)
				this->originindex->push_back((voxeltype)(coord[i] / increment[i]));//convert voxel center coordinate back to original voxel index
			for (; i<MaxDimension; i++)
				this->originindex->push_back(0);
		}
		//loop end
		coord[0] += increment[0];
		for (i = 0; i<this->VoxelDimension - 1; i++)
			if (coord[i]>finish)
			{
				coord[i] = initial[i];
				coord[i + 1] += increment[i+1];
			}
		if (coord[this->VoxelDimension - 1]<finish) goto start;

		//function end
		return;
	}
	void GetOriginalList(const CellList & Configuration, Shape * pShape)
	{
		short i;
		double increment[MaxDimension];
		for (i = 0; i<this->ParticleDimension; i++)
			increment[i] = 2.0*this->halfsize;
		for (; i<this->VoxelDimension; i++)
			increment[i] = 2.0*this->orienthalfsize;
		double finish = ::SystemSize;
		double initial[MaxDimension];
		for (i = 0; i < this->ParticleDimension; i++)
			initial[i] = this->halfsize;
		for (; i < this->VoxelDimension; i++)
			initial[i] = this->orienthalfsize;

		std::vector<double> coord0s;
		for (double temp = initial[0]; temp < finish; temp += increment[0])
			coord0s.push_back(temp);

#pragma omp parallel 
		{
			std::vector<voxeltype> myOriginindex;
			CellList myList(Configuration);
#pragma omp for schedule(guided)
			for (long ii = 0; ii < coord0s.size(); ii++)
			{
				double coord[MaxDimension];
				coord[0] = coord0s[ii];
				short i;
				//initialize
				for (i = 1; i < this->VoxelDimension; i++)
					coord[i] = initial[i];
				//main loop body
			start:
				if (myList.CheckVoxelOverlap(coord, pShape, this->halfsize, this->orienthalfsize) == false)
				{
					for (i = 0; i < this->VoxelDimension; i++)
						myOriginindex.push_back((voxeltype)(coord[i] / increment[i]));//convert voxel center coordinate back to original voxel index
					for (; i < MaxDimension; i++)
						myOriginindex.push_back(0);
				}
				//loop end
				coord[1] += increment[1];
				for (i = 1; i<this->VoxelDimension - 1; i++)
					if (coord[i]>finish)
					{
						coord[i] = initial[i];
						coord[i + 1] += increment[i+1];
					}
				if (coord[this->VoxelDimension - 1] < finish) goto start;
			}
			if (myOriginindex.size() > 0)
			{
#pragma omp critical
				{
					auto size = this->originindex->size();
					this->originindex->resize(size + myOriginindex.size());
					std::memcpy(&(*this->originindex)[size], &myOriginindex[0], sizeof(voxeltype)*myOriginindex.size());
				}
			}
		}

		//function end
		return;
	}
	unsigned long NbrVoxel(void) const
	{
		return this->originindex->size()/MaxDimension;
	}
	void GetCenter(double * result, size_t nbr) const//get the center coordinate of the nbr-th voxel
	{
		std::vector<voxeltype>::iterator itero=this->originindex->begin()+MaxDimension*nbr;
		std::vector<unsigned short>::iterator iters=this->subindex->begin()+this->level*nbr;
		double nowsize=this->originhalfsize;
		double orientNowsize=this->orientoriginhalfsize;
		int j;
		for(j=0; j<this->ParticleDimension; j++)
		{
			result[j]=(2*itero[j]+1)*nowsize;//result is center of original voxel
			assert(result[j]<=1+1e-4);
			assert(result[j]>=-1e-4);
		}
		for (; j<this->VoxelDimension; j++)
		{
			result[j] = (2 * itero[j] + 1)*orientNowsize;//result is center of original voxel
			assert(result[j] <= 1 + 1e-4);
			assert(result[j] >= -1e-4);
		}
		for (int i = 0; i<this->level; i++)
		{
			signed char * iterc=this->splitvoxellist+iters[i]*MaxDimension;
			nowsize/=2;
			orientNowsize/=2;
			for (j = 0; j<this->ParticleDimension; j++)
			{
				result[j] += nowsize*iterc[j];
				assert(result[j] <= 1 + 1e-4);
				assert(result[j] >= -1e-4);
			}
			for (; j<this->VoxelDimension; j++)
			{
				result[j] += orientNowsize*iterc[j];
				assert(result[j] <= 1 + 1e-4);
				assert(result[j] >= -1e-4);
			}
		}//result is the center of subvoxel
	}
	void GetRandomCoordinate(double * result, RandomGenerator & gen) const//generate a random coordinate which is inside a random voxel, write it into result
	{
		size_t nbr;
		do
		{
			nbr=(size_t)(gen.RandomDouble()*this->NbrVoxel());
		}
		while(nbr>=this->NbrVoxel());
		this->GetCenter(result, nbr);
		int j;
		for (j = 0; j<this->ParticleDimension; j++)
		{
			double rand2 = 2 * gen.RandomDouble() - 1;
			result[j] += rand2*this->halfsize;
			assert(result[j] <= 1 + 1e-4);
			assert(result[j] >= -1e-4);
		}
		for (; j<this->VoxelDimension; j++)
		{
			double rand2 = 2 * gen.RandomDouble() - 1;
			result[j] += rand2*this->orienthalfsize;
			assert(result[j] <= 1 + 1e-4);
			assert(result[j] >= -1e-4);
		}
	}
	bool CheckOverlapDeep(const CellList & Configuration, Shape * pShape, const double * Center, double QuaterSize, double orientQuaterSize, unsigned int checklevel)
	{
		if(Configuration.CheckVoxelOverlap(Center, pShape, 2*QuaterSize, 2*orientQuaterSize)==true) return true;
		if(checklevel==0)
			return false;
		else
		{
			double newCenter[::MaxDimension];
			for(unsigned short i=0; i<this->nbrsplitvoxel; i++)
			{
				int j;
				for(j=0; j< this->ParticleDimension; j++)
					newCenter[j]=Center[j]+QuaterSize*this->splitvoxellist[i*::MaxDimension+j];
				for(; j< this->VoxelDimension; j++)
					newCenter[j]=Center[j]+orientQuaterSize*this->splitvoxellist[i*::MaxDimension+j];
				if(Configuration.CheckOverlap(newCenter, pShape)==false)
					return false;
				if(this->CheckOverlapDeep(Configuration, pShape, newCenter, QuaterSize/2, orientQuaterSize/2, checklevel-1)==false)
					return false;
			}
			return true;
		}
	}
	void SplitVoxelSerial(const CellList & Configuration, Shape * pShape, int checklevel=0)
	{
		std::vector<voxeltype> * neworigin = new std::vector<voxeltype>;
		std::vector<unsigned short> * newsub = new std::vector<unsigned short>;
		neworigin->reserve(this->NbrVoxel()*MaxDimension);
		newsub->reserve(this->NbrVoxel()*(this->level+1));
		unsigned long nbr=this->NbrVoxel();

		for(signed long i=0; i<nbr; i++)
		{
			//locate this voxel
			std::vector<voxeltype>::iterator itero=this->originindex->begin()+MaxDimension*i;
			std::vector<unsigned short>::iterator iters=this->subindex->begin()+this->level*i;
			double center[MaxDimension], ncenter[MaxDimension];
			this->GetCenter(center, i);
			if(Configuration.CheckVoxelOverlap(center, pShape, this->halfsize, this->orienthalfsize)==true)
				continue;
			for(unsigned short j=0; j<this->nbrsplitvoxel; j++)
			{
				signed char * now=this->splitvoxellist+j*MaxDimension;
				//calculate new center
				int k;
				for (k = 0; k<this->ParticleDimension; k++)
					ncenter[k] = center[k] + now[k] * this->halfsize*0.5;
				for (; k<this->VoxelDimension; k++)
					ncenter[k] = center[k] + now[k] * this->orienthalfsize*0.5;
				if (this->CheckOverlapDeep(Configuration, pShape, ncenter, this->halfsize / 4, this->orienthalfsize / 4, checklevel) == false)
				{
					//add new voxel

					int k;
					for(k=0; k<this->VoxelDimension; k++)
						neworigin->push_back(itero[k]);
					for(; k<MaxDimension; k++)
						neworigin->push_back(0);
					for(k=0; k<this->level; k++)
						newsub->push_back(iters[k]);

					newsub->push_back(j);

				}
			}
		}

		delete this->originindex;
		delete this->subindex;
		this->originindex=neworigin;
		this->subindex=newsub;
		this->halfsize/=2;
		this->orienthalfsize /= 2;
		this->level++;

	}
	void SplitVoxel(const CellList & Configuration, Shape * pShape, int checklevel=0)
	{
		std::vector<voxeltype> * neworigin = new std::vector<voxeltype>;
		std::vector<unsigned short> * newsub = new std::vector<unsigned short>;
		neworigin->reserve(this->NbrVoxel()*MaxDimension);
		newsub->reserve(this->NbrVoxel()*(this->level+1));
		unsigned long nbr=this->NbrVoxel();

#pragma omp parallel default(shared)
		{
			CellList myList(Configuration);
			std::vector<voxeltype> myneworigin;
			std::vector<unsigned short> mynewsub;
#pragma omp for schedule(guided)
			for(signed long i=0; i<nbr; i++)
			{
				//locate this voxel
				std::vector<voxeltype>::iterator itero=this->originindex->begin()+MaxDimension*i;
				std::vector<unsigned short>::iterator iters=this->subindex->begin()+this->level*i;
				double center[MaxDimension], ncenter[MaxDimension];
				this->GetCenter(center, i);
				if (myList.CheckVoxelOverlap(center, pShape, this->halfsize, this->orienthalfsize) == true)
					continue;
				for(unsigned short j=0; j<this->nbrsplitvoxel; j++)
				{
					signed char * now=this->splitvoxellist+j*MaxDimension;
					//calculate new center
					int k;
					for (k = 0; k<this->ParticleDimension; k++)
						ncenter[k] = center[k] + now[k] * this->halfsize*0.5;
					for (; k<this->VoxelDimension; k++)
						ncenter[k] = center[k] + now[k] * this->orienthalfsize*0.5;
					if (this->CheckOverlapDeep(myList, pShape, ncenter, this->halfsize / 4, this->orienthalfsize/4, checklevel) == false)
					{
						//add new voxel
						for(k=0; k<this->VoxelDimension; k++)
							myneworigin.push_back(itero[k]);
						for(; k<MaxDimension; k++)
							myneworigin.push_back(0);
						for(k=0; k<this->level; k++)
							mynewsub.push_back(iters[k]);
							
						mynewsub.push_back(j);
					}
				}
			}
			if (myneworigin.size() > 0)
			{
#pragma omp critical
				{
					auto size = neworigin->size();
					neworigin->resize(size + myneworigin.size());
					std::memcpy(&(*neworigin)[size], &myneworigin[0], sizeof(voxeltype)*myneworigin.size());
					auto size2 = newsub->size();
					newsub->resize(size2 + mynewsub.size());
					std::memcpy(&(*newsub)[size2], &mynewsub[0], sizeof(unsigned short)*mynewsub.size());
				}
			}
		}

		delete this->originindex;
		delete this->subindex;
		this->originindex=neworigin;
		this->subindex=newsub;
		this->halfsize/=2;
		this->orienthalfsize /= 2;
		this->level++;
	}

	void Save(std::ostream & saveto) const
	{
		saveto.write((char*)(&this->level), sizeof(this->level));
		saveto.write((char*)(&this->halfsize), sizeof(this->halfsize));
		saveto.write((char*)(&this->originhalfsize), sizeof(this->originhalfsize));
		saveto.write((char*)(&this->orienthalfsize), sizeof(this->orienthalfsize));
		saveto.write((char*)(&this->orientoriginhalfsize), sizeof(this->orientoriginhalfsize));
		::SaveVectorData(*this->originindex, saveto);
		::SaveVectorData(*this->subindex, saveto);
	}
	void Load(std::istream & loadfrom)
	{
		loadfrom.read((char*)(&this->level), sizeof(this->level));
		loadfrom.read((char*)(&this->halfsize), sizeof(this->halfsize));
		loadfrom.read((char*)(&this->originhalfsize), sizeof(this->originhalfsize));
		loadfrom.read((char*)(&this->orienthalfsize), sizeof(this->orienthalfsize));
		loadfrom.read((char*)(&this->orientoriginhalfsize), sizeof(this->orientoriginhalfsize));
		::LoadVectorData(*this->originindex, loadfrom);
		::LoadVectorData(*this->subindex, loadfrom);
	}

	double GetVoxelSize() const
	{
		return this->halfsize*2;
	}

	void DistillVoxelList(const CellList & Configuration, Shape * pShape, int checklevel = 0)
	{
		std::vector<voxeltype> * neworigin = new std::vector<voxeltype>;
		std::vector<unsigned short> * newsub = new std::vector<unsigned short>;
		neworigin->reserve(this->NbrVoxel()*MaxDimension);
		newsub->reserve(this->NbrVoxel()*(this->level + 1));

		signed long max = this->NbrVoxel();
		progress_display pd(max);
#pragma omp parallel default(shared)
		{
			CellList myList(Configuration);
			std::vector<voxeltype> myneworigin;
			std::vector<unsigned short> mynewsub;
#pragma omp for schedule(guided)
			for (signed long i = 0; i < max; i++)
			{
#pragma omp critical(displayProgress)
				{
					pd++;
				}
				double center[MaxDimension];
				this->GetCenter(center, i);
				if (this->CheckOverlapDeep(myList, pShape, center, this->halfsize / 2, this->orienthalfsize / 2, checklevel))
					continue;
				for (int j = 0; j < ::MaxDimension; j++)
					myneworigin.push_back((*this->originindex)[i*::MaxDimension + j]);
				for (int j = 0; j < this->level; j++)
					mynewsub.push_back((*this->subindex)[i*this->level + j]);
			}
			if (myneworigin.size() > 0)
			{
#pragma omp critical
				{
					auto size = neworigin->size();
					neworigin->resize(size + myneworigin.size());
					std::memcpy(&(*neworigin)[size], &myneworigin[0], sizeof(voxeltype)*myneworigin.size());
					auto size2 = newsub->size();
					newsub->resize(size2 + mynewsub.size());
					std::memcpy(&(*newsub)[size2], &mynewsub[0], sizeof(unsigned short)*mynewsub.size());
				}
			}
		}
		delete this->originindex;
		delete this->subindex;
		this->originindex = neworigin;
		this->subindex = newsub;
	}

};


class Configuration
{
public:
	CellList * particles;
	VoxelList * voxels;
	RandomGenerator gen;
	double particlevolume;//the volume of a particle
	Shape * pShape;
	int ParticleDimension, VoxelDimension;

	Configuration(Shape * pShape)//NbrSphere: expected saturation number of particle
		: pShape(pShape)
	{
		//calculate diameter
		this->particlevolume = pShape->GetVolume();
		this->ParticleDimension = pShape->GetDimension();
		this->VoxelDimension = pShape->GetVoxelDimension();
		double InteractionRadius = pShape->GetInteractionRadius();
		unsigned int CellListCapacity = std::ceil(1.0 / this->particlevolume);

		particles=new CellList(pShape->GetDimension(), (CellListCapacity), InteractionRadius);
		voxels=nullptr;
	}
	~Configuration()
	{
		if(this->particles!=nullptr) delete this->particles;
		if(this->voxels!=nullptr) delete this->voxels;
	}

	void RSA_I(size_t NbrInsertedLimit, size_t TrialTime)//try to insert particle, if less than NbrInsertedLimit particles inserted in TrialTime trials, stop
	{
		size_t nbrinserted=10000;
		double temp[MaxDimension];
		while(nbrinserted>NbrInsertedLimit)
		{
			nbrinserted=0;
			for(size_t i=0; i<TrialTime; i++)
			{
				for(int j=0; j<this->VoxelDimension; j++)
					temp[j]=this->gen.RandomDouble();
				if(this->particles->CheckOverlap(temp, pShape)==false)
				{
					this->particles->Insert(temp);
					nbrinserted++;
				}
			}
		}
	}
	void RSA_II_Serial( size_t NbrInsertedLimit, size_t TrialTime)
	{
		if(this->voxels->NbrVoxel()==0) return;
		size_t nbrinserted=10000;
		double temp[MaxDimension];
		while(nbrinserted>NbrInsertedLimit)
		{
			nbrinserted=0;
			for(size_t i=0; i<TrialTime; i++)
			{
				this->voxels->GetRandomCoordinate(temp, this->gen);
				if(this->particles->CheckOverlap(temp, pShape)==false)
				{
					this->particles->Insert(temp);
					nbrinserted++;
				}
			}
		}
	}	
	struct RSA_II_TempStruct
	{
		double x[MaxDimension];
	};
	class RSA_II_TempClass
	{
	public:
		RandomGenerator * g;
		RSA_II_TempClass(RandomGenerator & gen)
		{
			this->g= & gen;
		}
		ptrdiff_t operator() (ptrdiff_t max)
		{
			return static_cast<ptrdiff_t>(this->g->RandomDouble()*max);
		}
	};
	void RSA_II_Parallel( size_t NbrInsertedLimit, size_t TrialTime)
	{
		if(this->voxels->NbrVoxel()==0) return;
		size_t nbrinserted=10000;
		std::vector<RSA_II_TempStruct> results;

		while(nbrinserted>NbrInsertedLimit)
		{
			nbrinserted=0;
#pragma omp parallel default(shared)
{
			CellList myList(*this->particles);
			RandomGenerator mygen((int)(this->gen.RandomDouble()*10000)+omp_get_thread_num());
			std::vector<RSA_II_TempStruct> myresults;
#pragma omp for
			for(signed long i=0; i<TrialTime; i++)
			{
				double temp[MaxDimension];
				this->voxels->GetRandomCoordinate(temp, mygen);
				if(myList.CheckOverlap(temp, pShape)==false)
				{
					myresults.push_back(RSA_II_TempStruct());
					for(int i=0; i< ::MaxDimension; i++)
						myresults.back().x[i]=temp[i];
					myList.Insert(temp);
				}
			}

#pragma omp critical
			{
				while(myresults.size()!=0)
				{
					results.push_back(myresults.back());
					myresults.pop_back();
				}
			}
}
			RSA_II_TempClass rtemp(this->gen);
			std::shuffle(results.begin(), results.end(), rtemp);
			for(auto iter=results.begin(); iter!=results.end(); iter++)
			{
				double * temp = iter->x;
				if(this->particles->CheckOverlap(temp, pShape)==false)
				{
					this->particles->Insert(temp);
					nbrinserted++;
				}
			}

			//std::cout<<"size of temp rsults:"<<results.size()<<'\n';

			results.clear();
		}
	}	
	void RSA_II(size_t NbrInsertedLimit, size_t TrialTime, bool parallel)
	{
		if(parallel)
			this->RSA_II_Parallel(NbrInsertedLimit, TrialTime);
		else
			this->RSA_II_Serial(NbrInsertedLimit, TrialTime);
	}
	double Density(void) const
	{
		return this->particles->nbrparticle * this->particlevolume;
	}
	void GetVoxelList(bool parallel, long InitialVoxelRank=0)//get original voxel list
	{
		if (InitialVoxelRank == 0)
		{
			//use default
			double radius = std::pow(this->particlevolume / UnitSphereVolume(this->ParticleDimension), 1.0 / this->ParticleDimension);
			InitialVoxelRank = std::ceil(2.0/radius);
		}
		this->voxels=new VoxelList(this->ParticleDimension, this->VoxelDimension, InitialVoxelRank, 24);
		if(parallel)
			this->voxels->GetOriginalList(* this->particles, pShape);
		else
			this->voxels->GetOriginalListSerial(* this->particles, pShape);
	}
	void SplitVoxel(bool parallel, int checklevel=0)
	{
		if(parallel)
			this->voxels->SplitVoxel(*this->particles, pShape, checklevel);
		else
			this->voxels->SplitVoxelSerial(*this->particles, pShape, checklevel);
	}

	void OutputSpheres(std::fstream & outfile)//output the coordinate of all particles
	{
		for(int i=1; i<this->particles->nbrparticle+1; i++)
		{
			for(int j=0; j<this->VoxelDimension; j++)
				outfile<<this->particles->coordinate[i*MaxDimension+j]<<'\t';
			outfile<<'\n';
		}
	}

	void Save(std::ostream & saveto) const
	{
		bool voxelsgenerated=(this->voxels!=nullptr);
		saveto.write((char*)(&voxelsgenerated), sizeof(voxelsgenerated));
		saveto.write((char*)(&this->ParticleDimension), sizeof(this->ParticleDimension));
		saveto.write((char*)(&this->VoxelDimension), sizeof(this->VoxelDimension));
		saveto.write((char*)(&this->particlevolume), sizeof(this->particlevolume));
		if(voxelsgenerated)
			this->voxels->Save(saveto);
		this->particles->Save(saveto);
	}
	void Load(std::istream & loadfrom)
	{
		bool voxelsgenerated;
		loadfrom.read((char*)(&voxelsgenerated), sizeof(voxelsgenerated));
		loadfrom.read((char*)(&this->ParticleDimension), sizeof(this->ParticleDimension));
		loadfrom.read((char*)(&this->VoxelDimension), sizeof(this->VoxelDimension));
		loadfrom.read((char*)(&this->particlevolume), sizeof(this->particlevolume));
		if(voxelsgenerated)
		{
			if(this->voxels!=nullptr)
				this->voxels->Load(loadfrom);
			else
			{
				this->voxels=new VoxelList(this->ParticleDimension, this->VoxelDimension, 10/*anything should work here*/, 10);
				this->voxels->Load(loadfrom);
			}
		}
		this->particles->Load(loadfrom);
	}
};



class Sphere : public Shape
{
private:
	int dimension;
	double radius;
	double Diameter2;
public:
	Sphere()
	{
		std::cout << "Sphere dimension=";
		std::cin >> dimension;
		std::cout << "Sphere radius=";
		std::cin >> radius;
		this->Diameter2 = 4 * radius*radius;
	}
	Sphere(int dimension, double radius) : dimension(dimension), radius(radius)
	{
		this->Diameter2 = 4 * radius*radius;
	}
	virtual int GetDimension()
	{
		return this->dimension;
	}
	virtual int GetVoxelDimension()
	{
		return this->dimension;
	}
	virtual double GetVolume()
	{
		return UnitSphereVolume(dimension)*std::pow(radius, dimension);
	}
	virtual bool Overlap(const double * adjParticle, const double * Coordinate)
	{
		double distance2 = 0;
		for (int j = 0; j < this->dimension; j++)
		{
			double temp = adjParticle[j] - Coordinate[j];
			distance2 += temp*temp;
		}
		if (distance2 < Diameter2)
			return true;
		else
			return false;
	}
	virtual bool VoxelOverlap(const double * adjParticle, const double * Coordinate, double HalfSize, double orientHalfSize)
	{
		double distance2 = 0;
		for (int j = 0; j<this->dimension; j++)
		{
			double temp = adjParticle[j] - Coordinate[j];
			if (temp>0)
				temp += HalfSize;
			else
				temp -= HalfSize;
			distance2 += temp*temp;
		}
		if (distance2<Diameter2)
			return true;
		else
			return false;
	}
	virtual double GetInteractionRadius()
	{
		return 2 * this->radius;
	}
};


#include <set>
class TwoDPolygon : public Shape
{
private:
	//polar coordinates of all vertices
	//assume vertex 0 is linked to vertex 1, vertex 1 is linked to vertex 2, vertex 2 is linked to vertex 3, etc.
	//assume vertex (VertexR.size()-1) is linked to vertex 0
	std::vector<double> VertexR, VertexTheta;

	//two polygons must overlap if they are within this distance
	double OverlapMinDist2;

	//if the polygon have parallel sides, a special check after voxel splitting is necessary
	bool HasParallelSides;

	void FindOverlapMinDist()
	{
		if (VertexR.size() < 3)
		{
			OverlapMinDist2 = 0.0;
			return;
		}
		double result = VertexR[0];
		for (int j = 0; j < this->VertexR.size(); j++)
		{
			size_t jj = j + 1;
			if (jj == this->VertexR.size())
				jj = 0;
			double x3 = VertexR[j] * std::cos(VertexTheta[j]);
			double y3 = VertexR[j] * std::sin(VertexTheta[j]);
			double x4 = VertexR[jj] * std::cos(VertexTheta[jj]);
			double y4 = VertexR[jj] * std::sin(VertexTheta[jj]);
			double t = (x4 - x3) / (y4 - y3);
			double d = std::abs(t*y3 - x3) / std::sqrt(1.0 + t*t);
			result = std::min(result, d);
		}
		OverlapMinDist2 = 4.0*result*result;
	}

	void FindHasParallelSides()
	{
		HasParallelSides = false;
		if (VertexR.size() < 3)
		{
			return;
		}
		for (int i = 0; i<this->VertexR.size(); i++)
			for (int j = 0; j < this->VertexR.size(); j++)
			{
				if (i == j)
					continue;
				size_t ii = i + 1;
				if (ii == this->VertexR.size())
					ii = 0;
				size_t jj = j + 1;
				if (jj == this->VertexR.size())
					jj = 0;

				double x1 =  VertexR[i] * std::cos(VertexTheta[i] );
				double y1 =  VertexR[i] * std::sin(VertexTheta[i] );
				double x2 =  VertexR[ii] * std::cos(VertexTheta[ii] );
				double y2 =  VertexR[ii] * std::sin(VertexTheta[ii] );
				double x3 =  VertexR[j] * std::cos(VertexTheta[j] );
				double y3 =  VertexR[j] * std::sin(VertexTheta[j] );
				double x4 =  VertexR[jj] * std::cos(VertexTheta[jj] );
				double y4 =  VertexR[jj] * std::sin(VertexTheta[jj] );

				double direction1 = std::atan2(y2 - y1, x2 - x1);
				double direction2 = std::atan2(y4 - y3, x4 - x3);
				double dDirection = std::abs(direction1 - direction2);
				if (dDirection < 1e-10 || std::abs(dDirection - Pi) < 1e-10)
				{
					HasParallelSides = true;
					return;
				}
			}
	}

	//calculate the area of the triangle made from the origin, vertex i, and vertex (i+1)
	double TriangleArea(size_t i)
	{
		size_t j = i + 1;
		if (j == VertexR.size())
			j = 0;

		//a, b, and c are the side lengths of the triangle
		double a = VertexR[i], b = VertexR[j];
		double dxc = VertexR[i] * std::cos(VertexTheta[i]) - VertexR[j] * std::cos(VertexTheta[j]);
		double dyc = VertexR[i] * std::sin(VertexTheta[i]) - VertexR[j] * std::sin(VertexTheta[j]);
		double c = std::sqrt(dxc*dxc + dyc*dyc);

		double s = 0.5*(a + b + c);
		return std::sqrt(s*(s - a)*(s - b)*(s - c));
	}

	//test if line segment from point 1 to 2 intersects with line segment from point 3 to 4
	static bool LineLineIntersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
	{
		double o1 = (y2 - y1)*(x3 - x2) - (x2 - x1)*(y3 - y2);
		double o2 = (y2 - y1)*(x4 - x2) - (x2 - x1)*(y4 - y2);
		double o3 = (y4 - y3)*(x1 - x4) - (x4 - x3)*(y1 - y4);
		double o4 = (y4 - y3)*(x2 - x4) - (x4 - x3)*(y2 - y4);

		return (o1*o2 < 0.0) && (o3*o4 < 0.0);
	}
	//same as above, except that endpoints 3 and 4 comes from a line in a voxel, and thus carry an uncertainty
	static bool LineVoxelIntersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4, double dx, double dtheta, double l3, double l4)
	{
		double o1 = (y2 - y1)*(x3 - x2) - (x2 - x1)*(y3 - y2);
		double o2 = (y2 - y1)*(x4 - x2) - (x2 - x1)*(y4 - y2);
		double o3 = (y4 - y3)*(x1 - x4) - (x4 - x3)*(y1 - y4);
		double o4 = (y4 - y3)*(x2 - x4) - (x4 - x3)*(y2 - y4);

		if ((o1*o2 < 0.0) && (o3*o4 < 0.0))
		{
			//the voxel center is intersecting
			double d3 = dx + dtheta*l3;
			double d4 = dx + dtheta*l4;
			double do1 = d3*(std::abs(y2 - y1) + std::abs(x2 - x1));
			double do2 = d4*(std::abs(y2 - y1) + std::abs(x2 - x1));
			double do3 = (d4 + d3)*(std::abs(x1 - x4) + std::abs(y1 - y4) + 2 * d4) + d4*(std::abs(y4 - y3) + std::abs(x4 - x3));
			double do4 = (d4 + d3)*(std::abs(x2 - x4) + std::abs(y2 - y4) + 2 * d4) + d4*(std::abs(y4 - y3) + std::abs(x4 - x3));

			if (do1 < std::abs(o1) && do2 < std::abs(o2) && do3 < std::abs(o3) && do4 < std::abs(o4))
				return true;
			else
				return false;
		}
		else
			return false;
	}
public:
	TwoDPolygon()
	{
		std::cout << "Enter a scaling factor for all vertices:\n";
		double factor;
		std::cin >> factor;
		std::cout << "Enter the polar coordinates of each vertex of the 2D polygon. Adjacent vertices are connected.  Data should be separated by only spaces, and should be, r1, theta1, r2, theta2, ... . Enter 0 0 to finish\n";
		for (;;)
		{
			double r, theta;
			std::cin >> r;
			std::cin >> theta;
			if (r > 0.0)
			{
				this->VertexR.push_back(r*factor);
				this->VertexTheta.push_back(theta);
			}
			else
				break;
		}
		FindOverlapMinDist();
		FindHasParallelSides();
	}
	virtual int GetDimension()
	{
		return 2;
	}
	virtual int GetVoxelDimension()
	{
		return 3;
	}
	virtual double GetVolume()
	{
		double result = 0.0;
		for (int i = 0; i < this->VertexR.size(); i++)
			result += this->TriangleArea(i);
		return result;
	}
	virtual double GetInteractionRadius()
	{
		double maxR = 0.0;
		for (int i = 0; i < this->VertexR.size(); i++)
			maxR = std::max(maxR, this->VertexR[i]);
		return 2.0*maxR;
	}
	virtual bool Overlap(const double * adjParticle, const double * Coordinate)
	{
		//easy check
		double distance2 = 0;
		for (int j = 0; j < 2; j++)
		{
			double temp = adjParticle[j] - Coordinate[j];
			distance2 += temp*temp;
		}
		if (distance2 < OverlapMinDist2)
			return true;

		//complex check
		for (int i = 0; i<this->VertexR.size(); i++)
			for (int j = 0; j < this->VertexR.size(); j++)
			{
				size_t ii = i + 1;
				if (ii == this->VertexR.size())
					ii = 0;
				size_t jj = j + 1;
				if (jj == this->VertexR.size())
					jj = 0;

				double x1 = adjParticle[0] + VertexR[i] * std::cos(VertexTheta[i] + 2 * Pi*adjParticle[2]);
				double y1 = adjParticle[1] + VertexR[i] * std::sin(VertexTheta[i] + 2 * Pi*adjParticle[2]);
				double x2 = adjParticle[0] + VertexR[ii] * std::cos(VertexTheta[ii] + 2 * Pi*adjParticle[2]);
				double y2 = adjParticle[1] + VertexR[ii] * std::sin(VertexTheta[ii] + 2 * Pi*adjParticle[2]);
				double x3 = Coordinate[0] + VertexR[j] * std::cos(VertexTheta[j] + 2 * Pi*Coordinate[2]);
				double y3 = Coordinate[1] + VertexR[j] * std::sin(VertexTheta[j] + 2 * Pi*Coordinate[2]);
				double x4 = Coordinate[0] + VertexR[jj] * std::cos(VertexTheta[jj] + 2 * Pi*Coordinate[2]);
				double y4 = Coordinate[1] + VertexR[jj] * std::sin(VertexTheta[jj] + 2 * Pi*Coordinate[2]);
				if (LineLineIntersect(x1, y1, x2, y2, x3, y3, x4, y4))
					return true;
			}
		return false;
	}
	virtual bool VoxelOverlap(const double * adjParticle, const double * Coordinate, double HalfSize, double orientHalfSize)
	{
		//easy check
		double distance2 = 0;
		for (int j = 0; j<2; j++)
		{
			double temp = adjParticle[j] - Coordinate[j];
			if (temp>0)
				temp += HalfSize;
			else
				temp -= HalfSize;
			distance2 += temp*temp;
		}
		if (distance2 < OverlapMinDist2)
			return true;

		//complex check
		for (int i = 0; i < this->VertexR.size(); i++)
			for (int j = 0; j < this->VertexR.size(); j++)
			{
				size_t ii = i + 1;
				if (ii == this->VertexR.size())
					ii = 0;
				size_t jj = j + 1;
				if (jj == this->VertexR.size())
					jj = 0;

				double x1 = adjParticle[0] + VertexR[i] * std::cos(VertexTheta[i] + 2 * Pi*adjParticle[2]);
				double y1 = adjParticle[1] + VertexR[i] * std::sin(VertexTheta[i] + 2 * Pi*adjParticle[2]);
				double x2 = adjParticle[0] + VertexR[ii] * std::cos(VertexTheta[ii] + 2 * Pi*adjParticle[2]);
				double y2 = adjParticle[1] + VertexR[ii] * std::sin(VertexTheta[ii] + 2 * Pi*adjParticle[2]);
				double x3 = Coordinate[0] + VertexR[j] * std::cos(VertexTheta[j] + 2 * Pi*Coordinate[2]);
				double y3 = Coordinate[1] + VertexR[j] * std::sin(VertexTheta[j] + 2 * Pi*Coordinate[2]);
				double x4 = Coordinate[0] + VertexR[jj] * std::cos(VertexTheta[jj] + 2 * Pi*Coordinate[2]);
				double y4 = Coordinate[1] + VertexR[jj] * std::sin(VertexTheta[jj] + 2 * Pi*Coordinate[2]);
				if (LineVoxelIntersect(x1, y1, x2, y2, x3, y3, x4, y4, HalfSize, 2 * Pi*orientHalfSize, VertexR[j], VertexR[jj]))
					return true;
			}
		return false;
	}
	virtual void ProcessSplittedVoxels(const CellList & cellList, VoxelList & voxels, std::ostream & output, RandomGenerator & gen)
	{
		//debug temp
		//return;

		if (HasParallelSides && voxels.NbrVoxel()>1000)
		{
			//select 100 voxels to test
			std::set<size_t> selectedVoxels;
			while (selectedVoxels.size() < 100)
				selectedVoxels.insert(std::floor(gen.RandomDouble()*voxels.NbrVoxel()));

			int maxCountWithinDistance = 0;
			for (auto iterj = selectedVoxels.begin(); iterj != selectedVoxels.end(); ++iterj)
			{
				int j = *iterj;
				double center0[::MaxDimension];
				voxels.GetCenter(center0, j);

				//sample 100 voxels, if 50 of them are very close to each other, then the parallel sides problem is suspected
				int countWithinDistance = 0;
				for (auto iteri = selectedVoxels.begin(); iteri != selectedVoxels.end(); ++iteri)
				{
					int i = *iteri;
					double center1[::MaxDimension];
					voxels.GetCenter(center1, i);
					double dx = center0[0] - center1[0];
					dx -= std::floor(dx + 0.5);
					double dy = center0[1] - center1[1];
					dy -= std::floor(dy + 0.5);
					double dist2 = dx*dx + dy*dy;
					if (dist2 < OverlapMinDist2)
						countWithinDistance++;
				}
				maxCountWithinDistance = std::max(maxCountWithinDistance, countWithinDistance);
			}

			if (maxCountWithinDistance>50)
			{
				//parallel sides problem suspected
				output << "at time" << std::time(nullptr) - start << ", " << "current voxel:" << voxels.NbrVoxel() << '\n';
				output << "within 100 sampled voxels, " << maxCountWithinDistance << " are very close to each other" << std::endl;
				output << "Parallel sides problem suspected. Running deep check.\n";
				output.flush();

				if(maxCountWithinDistance>99)
					voxels.DistillVoxelList(cellList, this, 12);
				else if(maxCountWithinDistance>95)
					voxels.DistillVoxelList(cellList, this, 6);
				else
					voxels.DistillVoxelList(cellList, this, 4);
			}
		}
	}
};



void SaveConfig(Configuration & con1, std::ostream & output)
{
			output<<"at time"<<std::time(nullptr)-start<<", Save Start";
			std::fstream ofile("System.dat", std::fstream::out | std::fstream::binary);
			con1.Save(ofile);
			std::fstream flag("continue.txt", std::fstream::out);
			flag<<"Needs to be Continued\n";
			ofile.close();
			flag.close();
			output<<"at time"<<std::time(nullptr)-start<<", Save Complete";
			exit(0);
}

double GenerateConfiguration(Shape * pShape, short nbrinsertedlimit, unsigned long trial1, unsigned long trial2, int seed, bool tryAgain, std::ostream & output, double * volumeratio=nullptr, bool Verbose=false)
{
	unsigned long OriginalTrial2=trial2;
	bool parallel=false;
	bool printconfig=true;
	if(nbrinsertedlimit<0)//use nbrinsertedlimit<0 to indicate parallel
	{
		nbrinsertedlimit*=-1;
		parallel=true;
	}

	Configuration con1(pShape);
	con1.gen.seed(seed);
	
	std::fstream ifile("System.dat", std::fstream::in | std::fstream::binary);
	if(ifile.good())
	{
		output<<"at time"<<std::time(nullptr)-start<<", Load Start";
		output.flush();
		con1.Load(ifile);
		output<<"at time"<<std::time(nullptr)-start<<", Load Complete";
		output.flush();
	}
	else
	{

		con1.RSA_I(nbrinsertedlimit, trial1);
		output<<"at time"<<std::time(nullptr)-start<<", "<<con1.particles->nbrparticle<<"particles inserted\n";
		output.flush();
		con1.GetVoxelList(parallel);
		if(Verbose)
		{
			output<<"voxel size:"<<con1.voxels->GetVoxelSize()<<" \tvoxel centers are:\n";
			for(size_t i=0; i<con1.voxels->NbrVoxel(); i++)
			{
				double temp[MaxDimension];
				con1.voxels->GetCenter(temp, i);
				for(size_t j=0; j<pShape->GetVoxelDimension(); j++)
					output<<temp[j]<<", \t";
				output<<",\n";
			}
			output.flush();
		}

	}
	ifile.close();
	int SplitCount=0;
	while(con1.voxels->NbrVoxel()!=0)
	{
		if(std::time(nullptr)-start>TimeLimit)
			SaveConfig(con1, output);
		output<<"at time"<<std::time(nullptr)-start<<", "<<"current voxel:"<<con1.voxels->NbrVoxel()<<'\n';
		output.flush();
		con1.RSA_II(nbrinsertedlimit, trial2, parallel);
		if(std::time(nullptr)-start>TimeLimit)
			SaveConfig(con1, output);
		output<<"at time"<<std::time(nullptr)-start<<", "<<con1.particles->nbrparticle<<"particles inserted\n";
		output.flush();

		unsigned long temp=con1.voxels->NbrVoxel();
		con1.SplitVoxel(parallel, 0);

		pShape->ProcessSplittedVoxels(*con1.particles, *con1.voxels, output, con1.gen);

		//double ratio=static_cast<double>(con1.voxels->NbrVoxel())/temp;
		//if(ratio<1)
		//	trial2*=ratio;
		//if(ratio>1.5)
		//	trial2=OriginalTrial2;
		SplitCount++;


		//debug temp
		//if (con1.voxels->NbrVoxel() > 5e7)
		//{
		//	//select 10 voxels to test
		//	RandomGenerator gen(23456);
		//	std::set<size_t> selectedVoxels;
		//	while (selectedVoxels.size() < 10)
		//		selectedVoxels.insert(std::floor(gen.RandomDouble()*con1.voxels->NbrVoxel()));

		//	std::cout << "Too many voxels, debugging\n";
		//	std::fstream outfile("SphereCoordinate.txt", std::ios::out);
		//	outfile.precision(17);
		//	con1.OutputSpheres(outfile);
		//	std::fstream outfile2("SelectedVoxels.txt", std::ios::out);
		//	outfile2.precision(17);

		//	for (auto iterj = selectedVoxels.begin(); iterj != selectedVoxels.end(); ++iterj)
		//	{
		//		int i = *iterj;

		//		double temp[MaxDimension];
		//		con1.voxels->GetCenter(temp, i);
		//		for (size_t j = 0; j<pShape->GetVoxelDimension(); j++)
		//			outfile2 << temp[j] << " \t";
		//		outfile2 << "\n";
		//	}
		//}

		if(Verbose)
		{
			output<<"voxel size:"<<con1.voxels->GetVoxelSize()<<" \tvoxel centers are:\n";
			for(size_t i=0; i<con1.voxels->NbrVoxel(); i++)
			{
				double temp[MaxDimension];
				con1.voxels->GetCenter(temp, i);
				for (size_t j = 0; j<pShape->GetVoxelDimension(); j++)
					output<<temp[j]<<", \t";
				output<<",\n";
			}
			output.flush();
		}
	}
	if(tryAgain)
	{
		delete con1.voxels;
		con1.GetVoxelList(parallel);
		while(con1.voxels->NbrVoxel()!=0)
		{
			output<<"at time"<<std::time(nullptr)-start<<", "<<"current voxel:"<<con1.voxels->NbrVoxel()<<'\n';
			output.flush();
			con1.SplitVoxel(parallel, 1);
		}
	}


	output<<"at time"<<std::time(nullptr)-start<<", configuration saturated\n";
	output<<"Volume of a Particle:"<<con1.particlevolume<<"\nDensity:"<<con1.Density()<<'\n';
	output.flush();
	if(volumeratio!=nullptr)
		*volumeratio=con1.particlevolume;

	if(printconfig)
	{
		std::fstream outfile("SphereCoordinate.txt", std::ios::out);
		outfile.precision(17);
		con1.OutputSpheres(outfile);
	}

	return con1.Density();
}


int cli()
{
	unsigned long nbrparticle, trial1, trial2;
	short nbrinsertedlimit;
	double voxelratio;
	int seed;

	std::cout<<"Enter the following\n\n1.InsertLimit 2.Trial1 and 3.Trial2\n(Remark: If less than InsertLimit particles inserted in Trial1 trials in 1st step, 1st step finish. The same for 2nd step.)\n4.Random Seed\n5.TimeLimit\n";
	std::cin>>nbrinsertedlimit;
	std::cin>>trial1;
	std::cin>>trial2;
	std::cin>>seed;
	std::cin>>TimeLimit;

	std::unique_ptr<Shape> pShape;
	std::cout << "Choose Shape\n\n1.Sphere\n2.2D Polygon\nchoice=";
	int ShapeType;
	std::cin >> ShapeType;
	if (ShapeType == 1)
	{
		pShape = std::unique_ptr<Shape>(new Sphere);
	}
	else if (ShapeType == 2)
	{
		pShape = std::unique_ptr<Shape>(new TwoDPolygon);
	}
	else
	{
		std::cerr << "Unsupported shape!\n";
		return 1;
	}

	
	GenerateConfiguration(pShape.get(), nbrinsertedlimit, trial1, trial2, seed, false, std::cout);

	std::fstream flag("finish.txt", std::fstream::out);
	flag<<"Finished\n";

	return 0;
}

int main()
{
	start=std::time(nullptr);
	std::cout.precision(17);
	return cli();
}
