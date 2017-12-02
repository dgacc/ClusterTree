package Structures;

import java.util.ArrayList;
import java.util.Random;

import Files_InOut.ReadFiles;
import Operator.Evaluation;
import Operator.Selection;

public class Population {
	public static int populationLength = 100; // default length of population
	public ArrayList<Individual> population = new ArrayList<Individual>(populationLength);
    private Selection selection = new Selection();
    Individual individual = new Individual();
    Evaluation evaluation = new Evaluation();
    
	
	public Population() {
	};
	
	// get gene of an Individual
	public Individual getIndividual(int i) {
		return population.get(i);
	}
	
	public void addIndiv(Individual individual){
		population.add(individual);
	}
	public double[] populationFitness(){
		double[] popFitness = new double[populationLength];
		for( int i = 0; i < populationLength; i++){
			popFitness[i] = evaluation.evaluation(ReadFiles.weightMatrix, population.get(i).getGene(), ReadFiles.num_vertex, ReadFiles.root);
		}
		return popFitness;
	}
	public double[] populationFitness1(){
		double[] popFitness = new double[populationLength];
		for( int i = 0; i < populationLength; i++){
			popFitness[i] = evaluation.distanceEvaluate(ReadFiles.weightMatrix, population.get(i).getGene(), ReadFiles.num_vertex, ReadFiles.root);
		}
		return popFitness;
	}
	
	}
