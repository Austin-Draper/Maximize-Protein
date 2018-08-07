///////////////////////////////////////////////////////////////////////////////
// maxprotein.hh
//
// Compute the set of foods that maximizes protein, within a calorie budget,
// with the greedy method or exhaustive search.
//
///////////////////////////////////////////////////////////////////////////////


#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <queue>
#include <sstream>
#include <string>
#include <vector>
#include "timer.hh"
using namespace std;
// Simple structure for a single protein
struct Protein {
	Protein() {
		description = "";
		sequence = "";
	}
	Protein(std::string desc, std::string seq) {
		description = desc;
		sequence = seq;
	}
std::string		description;
std::string 	sequence;
};

// Alias for a vector of shared pointers to Protein objects.
typedef std::vector<std::shared_ptr<Protein>> ProteinVector;


// -------------------------------------------------------------------------
// Load all the proteins from a standard FASTA format file with one line
// per sequence (multi-line sequences are not allowed).
// Returns false on I/O error.
// -------------------------------------------------------------------------
bool load_proteins(ProteinVector & proteins, const std::string& path)
{
	//std::cout << "Loading proteins from [" << path << "]" << std::endl;
	proteins.clear();
	std::ifstream ifs(path.c_str());
	if (!ifs.is_open() || !ifs.good()) {
		std::cout << "Failed to open [" << path << "]" << std::endl;
		return false;
	}
	int proteinsLoaded = 0;
	bool have_description = false;
	std::shared_ptr<Protein> newProtein = nullptr;
	while (!ifs.eof()) {
		std::string lineBuffer;
		std::getline(ifs, lineBuffer);
		if (ifs.eof()) {
			break;
		}
		if (lineBuffer.size() == 0) {
			continue;
		}
		if (lineBuffer[0] == '>') {
			newProtein = std::shared_ptr<Protein>(new Protein);
			newProtein->description = lineBuffer.substr(1);
			have_description = true;
		}
		else if (have_description) {
			newProtein->sequence = lineBuffer;
			proteins.push_back(newProtein);
			proteinsLoaded++;
			have_description = false;
		}
	}

	ifs.close();
	//std::cout << "Loaded " << proteinsLoaded << " proteins from [" << path << "]" << std::endl;

	return true;
}


// -------------------------------------------------------------------------
int dynamicprogramming_longest_common_subsequence(const std::string & string1,
	const std::string & string2)
{
	int n = string1.size();
	int m = string2.size();
	vector<vector<int>> array2D(n+1, vector<int>(m+1, 0));
	
	if (n == 0 || m == 0)
		return 0;
	if (n == 1 || m == 1)
	{
		if (string1[0] == string2[0])
			return 1;
		return 0;
	}
	for (int i = 1; i <= n; i++)
	{
		for (int j = 1; j <= m; j++)
		{
			int up = array2D[i - 1][j];
			int left = array2D[i][j - 1];
			int diag = array2D[i - 1][j - 1];
			if (string1[i - 1] == string2[j - 1])
				diag = diag + 1;
			array2D[i][j] = max(up, left);
			array2D[i][j] = max(array2D[i][j], diag);
		}
	}
	
	return array2D[n][m];
}


// -------------------------------------------------------------------------
std::unique_ptr<std::vector<std::string>> generate_all_subsequences(const std::string & sequence)
{
	//std::unique_ptr<std::vector<std::string>> subsequences(nullptr);
	string subsequence = "";
	//R = {}
	int str_size = sequence.size();
	int n = pow(sequence.size(), 2);
	std::unique_ptr<std::vector<std::string>> subsequences(new vector<string>(sequence.size()));

	if (sequence.size() == 1)
	{
		//cout << "size ended up being 1\n";
		subsequence.append(sequence, 0);
		cout << subsequence;
		subsequences->push_back(subsequence);  
		//cout << "returning the size 1 subsequence\n";
		return subsequences;
	}
	for (int x = 0; x <= n - 1; x++)
	{
		subsequence = "";
		for (int j = 0; j <= str_size - 1; j++)
		{
			if (((x >> j) & 1) == 1)
			{
				//cout << "Sequence j is: " << sequence[j] << endl;
				subsequence = subsequence + sequence[j];
			}
			
		}
		//cout << "current subsequence is: " << subsequence << "\n";
		subsequences->push_back(subsequence);
		//subsequences->insert(subsequences->begin() + x, subsequence);
	}
	//cout << "returning subsequence" << endl;
	return subsequences;
					
}


// -------------------------------------------------------------------------
int exhaustive_longest_common_subsequence(const std::string & string1,
	const std::string & string2)
{
	std::unique_ptr<std::vector<std::string>> all_subseqs1(nullptr);
	std::unique_ptr<std::vector<std::string>> all_subseqs2(nullptr);
	all_subseqs1 = generate_all_subsequences(string1);
	all_subseqs2 = generate_all_subsequences(string2);
	string s1;
	string s2;
	int best_score = 0;
	int counter = 1;

	for (int i = 0; i < all_subseqs1->size(); i++)
	{
		for (int j = 0; j < all_subseqs2->size(); j++)
		{
			s1 = all_subseqs1->at(i);
			s2 = all_subseqs2->at(j);
			if (s1 == s2 && s1.length() > best_score)
			{
				best_score = s1.length();
			}
		}
	}

	return best_score;

}


// -------------------------------------------------------------------------
std::shared_ptr<Protein> exhaustive_best_match(ProteinVector & proteins, const std::string & string1)
{

	std::shared_ptr<Protein> best_protein(new Protein);
	int best_i = 0;
	int best_score = 0;
	string protein;
	string protein_desc;//delete
	string best_protein_str;
	string best_protein_desc;//delete
	int str_size = proteins.size();
	for (int i = 0; i <= str_size - 1; i++)
	{
		protein = proteins.at(i)->sequence;
		protein_desc = proteins.at(i)->description;//delete
		int score = exhaustive_longest_common_subsequence(protein, string1);
		if (score > best_score)
		{
			best_score = score;
			best_i = i;
			best_protein_str = protein;
			best_protein_desc = protein_desc;//delete
		}
	}
	best_protein->sequence = best_protein_str;
	cout << "\n" << best_protein_desc << "\n";//delete
	return best_protein;
}


// -------------------------------------------------------------------------
std::shared_ptr<Protein> dynamicprogramming_best_match(ProteinVector & proteins, const std::string & string1)
{
	std::shared_ptr<Protein> best_protein(new Protein);
	int best_i = 0;
	int best_score = 0;
	string protein;
	string protein_desc;//delete
	string best_protein_str;
	string best_protein_desc;//delete
	int str_size = proteins.size();
	cout << "str_size is: " << str_size << endl;
	for (int i = 0; i <= str_size - 1; i++)
	{
		protein = proteins.at(i)->sequence;
		protein_desc = proteins.at(i)->description;//delete
		int score = dynamicprogramming_longest_common_subsequence(protein, string1);

		if (score > best_score)
		{
			best_score = score;
			best_i = i;
			best_protein_str = protein;
			best_protein_desc = protein_desc;//delete
		}
	}
	best_protein->sequence = best_protein_str;
	cout << "\n" << best_protein_desc << "\n";//delete
	return best_protein;
}
