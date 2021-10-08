#include "Fossils.hh"

Fossils::Fossils()
{
}

Fossils::~Fossils()
{
	i_fos iF;
	iF = FossilRecord.begin();
	while (iF != FossilRecord.end())
	{
		if (*iF != NULL)
		{
			delete (*iF);
			iF = FossilRecord.erase(iF);
		}
		else	iF++;
	}
}

void Fossils::EraseFossil(unsigned long long fossilID)
{
	i_fos iF;
	iF = FossilRecord.begin();
	while (iF != FossilRecord.end())
	{
		if ((*iF)->fossil_id == fossilID)
		{
			iF = FossilRecord.erase(iF);
			return;
		}
		iF++;
	}
}

void Fossils::BuryFossil(Organelle* O)
{
	FossilRecord.push_back(O);
}

bool compare_fossils(const Organelle* first, const Organelle* second)
{
    if(first->fossil_id < second->fossil_id) return true;
    else if(first->fossil_id > second->fossil_id) return false;
		return true;	//Impossible!
}

void Fossils::SortFossils()
{
	FossilRecord.sort(compare_fossils);
}

void Fossils::ExhibitFossils()
{
	FILE* f;
	char OutputFile[800];
	i_fos iF;
	sprintf(OutputFile, "%s/ancestors/anctrace%08d.txt", folder.c_str(), Time);
	f=fopen(OutputFile, "w");

	if (f == NULL)	printf("Failed to open file for writing the ancestor trace.\n");
	// fprintf(f, "#id\t#anc_id\t#time_oa\t#genome\n");	//Don't print the header to save space, but at least you know now!

	iF = FossilRecord.begin();
	while(iF != FossilRecord.end())
	{
		if ((*iF)->Ancestor == NULL)
		{
			fprintf(f, "%llu\t%d\t%d\t%f\t%s\t%s\n", (*iF)->fossil_id, 0, (*iF)->time_of_appearance, (*iF)->nutrient_claim, (*iF)->G->ShowDefinition(false).c_str(), (*iF)->G->Show(NULL, false, true).c_str());
		}
		else
		{
			fprintf(f, "%llu\t%llu\t%d\t%f\t%s\t%s\n", (*iF)->fossil_id, ((*iF)->Ancestor)->fossil_id, (*iF)->time_of_appearance, (*iF)->nutrient_claim, (*iF)->G->ShowDefinition(false).c_str(), (*iF)->G->Show(NULL, false, true).c_str());
		}
		iF++;
	}
	fclose(f);
}
