package Cluster_Tree;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
//import java.util.ArrayList;
import java.util.Random;
import javax.swing.*;

import Files_InOut.ReadFiles;
import Operator.Crossover;
import Operator.Mutations;
import Operator.Selection;
import Structures.Individual;
import Structures.Population;
public class GA{
	private static double m_rate = 0.15;
	private static final double p_rate = 0.8; // crossover rate
	private static  Random r = new Random();  
	private static  int generation = 5000;  
	public static int populationLength = 100; // default length of population
	private static Crossover crossover = new Crossover();
	private static Mutations mutation = new Mutations();
	private static Selection selection = new Selection();

	
	public static void main(String[] args) {
		ReadFiles.clusterReadFiles("C:/Users/TrungTB/Desktop/test/5eil51.clt");
//	    ReadFiles.clusterReadFiles("C:/Users/TrungTB/Desktop/test/" +args[0]+ ".clt");
		int  num_vertex = ReadFiles.num_vertex;
	
//		JFrame gf = new JFrame();
//		gf.setVisible(true);
//		gf.setSize(800, 800);
//		gf.setTitle("the best Individual");
//		gf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//		gf.setVisible(true);
	
///*-----------------------------------------Build the new instance--------------------------------------------------------------*/
//		//optimal tree
//		 double[][] optimalTree = new double[num_vertex][num_vertex];
//		 // build the weight matrix from optimalTree
//		 double[][] weightMatrix1 = new double[num_vertex][num_vertex]; 
//		 Individual individual1 = new Individual();
//		 individual1.inintilizeIndividual();
//		 optimalTree = individual1.getGene();
//		 ArrayList<Cluster> clusters = ReadFiles.clusters;
//		 
//		 for(int i = 0; i < clusters.size(); i ++){
//				int numberClusterVertex = clusters.get(i).getCluster().size();
//				for( int j = 0; j < numberClusterVertex; j++ ){
//					for( int k = 0; k < numberClusterVertex; k++){
//						double value = 11 + r.nextInt(10);
//						weightMatrix1[clusters.get(i).getCluster().get(j)][clusters.get(i).getCluster().get(k)] = value;
//						weightMatrix1[clusters.get(i).getCluster().get(k)][clusters.get(i).getCluster().get(j)] = value;
//					}
//				}
//			}
//		 for( int i = 0; i < num_vertex; i++){
//			 for(int j = 0; j < num_vertex; j++){
//				if(optimalTree[i][j] > 0){
//				weightMatrix1[i][j] = optimalTree[i][j];
//				
//				}
//			 }
//		 }
//		 for( int i = 0; i < num_vertex; i++){
//			 for(int j = 0; j < num_vertex; j++){
//				if(weightMatrix1[i][j] == 0){
//					double value = 26 + r.nextInt(75);
//				    weightMatrix1[i][j] = value; 
//				}
//			 }
//		 }
//		 
//		 for( int i = 0; i < num_vertex; i ++){
//			 for( int j = 0; j < num_vertex; j++){
//				System.out.print(" " +weightMatrix1[i][j]);
//			 }
//			 System.out.println();
//		 }
//		 
//		 Evaluation evaluation = new Evaluation();
////		 Paint  p = new Paint();
////			p.weightMatrix = optimalTree;
////			p.fitness = evaluation.distanceEvaluate(weightMatrix1, optimalTree, num_vertex, ReadFiles.root);
////			p.num_vertex = num_vertex;
////		    gf.add(p);
////		    gf.setVisible(true);
//		
//		   // print files:
//		    CreateFile g = new CreateFile();
//		    g.cost =  evaluation.distanceEvaluate(weightMatrix1, optimalTree, num_vertex, ReadFiles.root);
//		    g.weightMatrix = weightMatrix1;
//	     	g.writeFile();
//	     	
//		    
///*--------------------------------------------------------END-------------------------------------------------------------*/
		  PrintWriter pw = null;
			try {
				pw = new PrintWriter(new FileWriter(new File(args[1]), true));
			} catch (IOException e) {
				e.printStackTrace();
			}
		 
		  
		  for(int loop = 0; loop < 7; loop ++){
			  
		    System.out.println(loop);
			// Initialize population
		Population population = new Population();
		for( int i = 0; i < populationLength; i++  ){
	        Individual individual = new Individual();
	        individual.inintilizeIndividual();
			population.addIndiv(individual);
		}
		
		System.out.print("---------------------------------------------> CLUSTER TREE <--------------------------------------------------\n ");
		System.out.println("Crossover rate : " + p_rate);
		System.out.println("Mutation  rate : " + m_rate);
		System.out.println("The number of generation : " +generation);
		
		for (int i = 0; i <  generation; i++) {
			Population subPop = new Population(); // store new individuals are generated
			double[] popFitness = population.populationFitness();
//			double[] popFitness1 = population.populationFitness1();
			
			
			for (int j = 0; j < Population.populationLength / 2; j++) {
				Individual father = new  Individual();
				Individual mother = new  Individual();
				
				//use tournament to  select Individual 
				father = population.getIndividual(selection.touranmentSelection(popFitness, 2, r));
				mother = population.getIndividual(selection.touranmentSelection(popFitness, 2, r));		
				
					
				Individual offspring1 =  new Individual();
				Individual offspring2 =  new Individual();
			    
				// new  random number  [0;1]
				double d = 0 + (1 - 0) * r.nextDouble();
				// if random number >  crossover rate,  add  the individual to population without crossover 
		     	if (d > p_rate) {
					subPop.addIndiv(father);
					subPop.addIndiv(mother);
			     } else {
			//	 else do crossover, generate two offsprings, then  do mutation
			    	offspring1.setGene(crossover.clusterCrossover(father.getGene(), mother.getGene(),
	   		    		       num_vertex, ReadFiles.clusters));
			    	offspring2.setGene(crossover.clusterCrossover(father.getGene(), mother.getGene(),
		    		       num_vertex, ReadFiles.clusters)); 
//					offspring1.setGene(crossover.ClusterBFSCrossover(father.getGene(), mother.getGene(),
//   		    		       num_vertex, ReadFiles.clusters));
//					
//					offspring2.setGene(crossover.ClusterBFSCrossover(father.getGene(), mother.getGene(),
//	   		    		       num_vertex, ReadFiles.clusters));
//					
					if( d < m_rate){
						 offspring1.setGene(mutation.mutationClusterTree(mother.getGene(),
								 num_vertex, ReadFiles.clusters));
						 offspring2.setGene(mutation.mutationClusterTree(father.getGene(), 
								 num_vertex, ReadFiles.clusters));
					}
					
					
				subPop.addIndiv(offspring1);			 
				subPop.addIndiv(offspring2); 	
				
			}
			}
			// get the index of the best Individual  of old Population
			double bestFiness = popFitness[0];
			int bestFinessIndex = 0;
			
			for(int t = 1; t < populationLength; t++ ){
			 if(popFitness[t] < bestFiness){
				 bestFiness = popFitness[t];
				  bestFinessIndex = t;
			 }
			
			}
			// paint graph	
//			Paint  p = new Paint();
//			p.weightMatrix = population.getIndividual(bestFinessIndex).getGene();
//			p.fitness = popFitness[bestFinessIndex];
//			p.fitness1 = popFitness1[bestFinessIndex];
//			p.num_vertex = num_vertex;
//		    gf.add(p);
//		    gf.setVisible(true);
//		    
//		    
			// add the best Individual  of old Population to new Population
			subPop.population.set(r.nextInt(Population.populationLength), 
					population.getIndividual(bestFinessIndex));	
					population = subPop;
		}
		
		double[] popFitness = population.populationFitness();
		double bestFiness = popFitness[0];
		int bestFinessIndex = 0;
		for(int t = 1; t < populationLength; t++ ){
			 if(popFitness[t] < bestFiness){
				 bestFiness = popFitness[t];
				  bestFinessIndex = t;
			 }
			}
//		
//		JFrame gf1 = new JFrame();
//		gf1.setVisible(true);
//		gf1.setSize(800, 800);
//		gf1.setTitle("the best ever");
//		gf1.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//		gf1.setVisible(true);
//		Paint  p = new Paint();
//		p.weightMatrix = population.getIndividual(bestFinessIndex).getGene();
//		p.fitness = popFitness[bestFinessIndex];
//		p.num_vertex = num_vertex;
//	    gf1.add(p);
//	    gf1.setVisible(true);
	    
	    // print to files
		 pw.println(popFitness[bestFinessIndex]);
		 
		 
	    
		
		// the best Individual at all
	   //	System.out.print("The best prufer number after " + generation + " generations: \n");
	}
		  pw.close();
		  }
}

