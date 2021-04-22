//
// Created by Rina Karnauch & Ofek Kaveh on 14/04/2021.
//


/**
 * main for the struct alignment
 * @param argc amount of args to the program
 * @param argv array of arguments as chars
 * @return none
 */
#include "Vector3.h"
#include "Atom.h"
#include "RigidTrans3.h"
#include "Matrix3.h"
#include "Molecule.h"
#include "PDB.h"
#include "Match.h"
#include "GeomHash.h"
#include <chrono>
#include <iostream>
#include "Triangle.h"

using namespace std;


vector<Triangle> makeCombinations(const Molecule<Atom> &molecule)
{
	vector<Triangle> triangles;
	for (int i = 0; i < molecule.size() - 3; i++)
	{

		Atom x = molecule[i];
		Atom y = molecule[i + 1];
		Atom z = molecule[i + 2];
		Triangle t = Triangle(x, y, z);
		triangles.push_back(t);
	}
	return triangles;
}

RigidTrans3 getTransformation(Triangle model, Triangle target)
{
	return target | model;
}

vector<RigidTrans3>
getAllTransformations(vector<Triangle> mulModel,
					  vector<Triangle> mulTarget)
{
	vector<RigidTrans3> transformations;
	for (int i = 0; i < mulModel.size(); i++)
	{
		for (int j = 0; j < mulTarget.size(); j++)
		{
			transformations.push_back(getTransformation(mulModel[i], mulTarget[j]));
		}
	}
	return transformations;

}

int main(int argc, char *argv[])
{
	// measure the run time
	auto start = std::chrono::system_clock::now();

	if (argc != 4)
	{
		std::cerr << "Usage: " << argv[0] << " target.pdb model.pdb dist_threshold"
				  << std::endl;
		exit(1);
	}

	// ********************Parameters********************
	float m_fDistThr = atof(argv[3]); // distance threshold on atoms in correspondence

	std::cout << "Distance threshold: " << m_fDistThr << std::endl;

	// read the two files into Molecule
	Molecule<Atom> molModel, molTarget, molModelFull;

	std::ifstream fileModel(argv[2]);
	std::ifstream fileModelFull(argv[2]);
	std::ifstream fileTarget(argv[1]);

	if (!fileModel)
	{
		std::cout << "File " << argv[1] << "does not exist." << std::endl;
		return 0;
	}
	if (!fileTarget)
	{
		std::cout << "File " << argv[2] << "does not exist." << std::endl;
		return 0;
	}

//
	int modelProteinSize = molModel.readPDBfile(fileModel, PDB::CAlphaSelector());
	int modelTargetSize = molTarget.readPDBfile(fileTarget, PDB::CAlphaSelector());
	molModelFull.readPDBfile(fileModelFull, PDB::AllSelector());


//	int modelProteinSize = molModel.readPDBfile(fileModel, PDB::AllSelector());
//	int modelTargetSize = molTarget.readPDBfile(fileTarget, PDB::AllSelector());

	// one is RNA and second is Protein or otherwise
	if (((modelProteinSize == 0) && (modelTargetSize > 0)) ||
		((modelProteinSize) > 0 && (modelTargetSize == 0)))
	{
		std::cout << "Alignment Not Possible between Protein and RNA" << std::endl;
		return 0;
	}
		// if no CA's found then its a RNA mol
	else if ((modelProteinSize == 0) && (modelTargetSize == 0))
	{
		modelProteinSize = molModel.readPDBfile(fileModel, PDB::PSelector());
		modelTargetSize = molTarget.readPDBfile(fileTarget, PDB::PSelector());
	}

	if ((modelProteinSize == 0) && (modelTargetSize == 0))
	{
		std::cout << "No Molecules Found." << std::endl;
		return 0;
	}

	// calculate center of mass
	Vector3 vectModelMass(0, 0, 0);
	Vector3 vectModelMassFull(0, 0, 0);
	for (unsigned int i = 0; i < molModel.size(); i++)
	{
		vectModelMass += molModel[i].position();
	}
	vectModelMass /= molModel.size();

	Vector3 vectTargetMass(0, 0, 0);
	Vector3 vectTargetMassFull(0, 0, 0);
	for (unsigned int i = 0; i < molTarget.size(); i++)
	{
		vectTargetMass += molTarget[i].position();
	}
	vectTargetMass /= molTarget.size();

	// transform the molecules to the center of the coordinate system
	molModel += (-vectModelMass);
	molTarget += (-vectTargetMass);

	// next we insert the target molecule into hash
	// this will help us to find atoms that are close faster
	GeomHash<Vector3, int> gHash(3,
								 m_fDistThr); // 3 is a dimension and m_fDistThr is the size of the hash cube
	for (unsigned int i = 0; i < molTarget.size(); i++)
	{
		gHash.insert(molTarget[i].position(),
					 i); // coordinate is the key to the hash, we store atom index
	}

	// now molecules are both from the same genre
	vector<Triangle> modelTriangles = makeCombinations(molModel);
	vector<Triangle> targetTriangles = makeCombinations(molTarget);

	vector<RigidTrans3> transformations = getAllTransformations(modelTriangles, targetTriangles);

	// now we try rotations and choose the best alignment from random rotations
	unsigned int iMaxSize = 0;
	RigidTrans3 rtransBest = transformations[0];
	float RMSD = 0.0;

	for (int it = 0; it < transformations.size(); it++)
	{

		// match is a class that stores the correspondence list, eg.
		// pairs of atoms, one from each molecule, that are matching
		RigidTrans3 trans = transformations[it];
		Match match = Match(trans);
		Matrix3 transformation_matrix = trans.rotation();

		// apply rotation on each atom in the model molecule and
		// add the pairs of atoms (one from target and one from model)
		// that are close enough to the match list
		for (unsigned int i = 0; i < molModel.size(); i++)
		{
			Vector3 mol_atom = transformation_matrix * molModel[i].position(); // rotate

			// find close target molecule atoms using the hash
			HashResult<int> result;
			gHash.query(mol_atom, m_fDistThr, result); // key is mol atom coordinate

			// check if the atoms in the result are inside the distance threshold
			// the hash is a cube shape, there can be atoms further that the threshold
			for (auto x = result.begin(); x != result.end(); x++)
			{
				float dist = mol_atom.dist(molTarget[*x].position());
				if (dist <= m_fDistThr)
				{
					float score = (1 / (1 + dist));
					match.add(*x, i, score, score);
				}
			}
			result.clear();
		}

		match.calculateBestFit(molTarget, molModel);

		if (iMaxSize < match.size())
		{
			iMaxSize = match.size();
			rtransBest = match.rigidTrans();
			RMSD = match.rmsd();
		}
	}

	std::cout << "Max Alignment Size: " << iMaxSize << std::endl;
	std::cout << "Achieved RMSD is:" << RMSD << std::endl;
	std::cout << "Rigid Trans: " << rtransBest << std::endl;

	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end - start;
	std::cout << "elapsed time: " << elapsed_seconds.count() << "s" << std::endl;


	molModelFull *=rtransBest;

	std::ofstream ofstream("transformed.pdb", std::ofstream::out);
	ofstream << molModelFull;
}


