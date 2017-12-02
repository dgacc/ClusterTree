package Cluster_Tree;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import javax.swing.*;

import Dislay.Paint;
import Files_InOut.ReadFiles;
import Operator.Crossover;
import Operator.Mutations;
import Operator.Selection;
import Structures.Individual;
import Structures.Population;
import random.MyRandom;
public class GA{
	private static double m_rate = 0.15;
	private static final double p_rate = 0.8; // crossover rate
	private static  int generation = 10000;  
	public static int populationLength = 100; // default length of population
	private static Crossover crossover = new Crossover();
	private static Mutations mutation = new Mutations();
	private static Selection selection = new Selection();


	
	public static void main(String[] args) {
		ReadFiles.clusterReadFiles("C:/Users/TrungTB/Desktop/test/10eil51.clt");
//	    ReadFiles.clusterReadFiles("C:/Users/TrungTB/Desktop/test/" +args[0]+ ".clt");
		int  num_vertex = ReadFiles.num_vertex;
////	
		JFrame gf = new JFrame();
		gf.setVisible(true);
		gf.setSize(800, 800);
		gf.setTitle("the best Individual");
		gf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		gf.setVisible(true);
//	

//		  PrintWriter pw = null;
//			try {
//				pw = new PrintWriter(new FileWriter(new File(args[1]), true));
//			} catch (IOException e) {
//				e.printStackTrace();
//			}
//		 
		  
		  for(int loop = 0; loop < 1; loop ++){
			MyRandom.setSeed(loop);  
		    System.out.println(loop);
			// Initialize population
		Population population = new Population();
		for( int i = 0; i < populationLength; i++  ){
	        Individual individual = new Individual();
	        individual.inintilizeIndividual();
			population.addIndiv(individual);
		}
		
		System.out.print("---------------------------------> CLUSTER TREE <--------------------------------------\n ");
		System.out.println("Crossover rate : " + p_rate);
		System.out.println("Mutation  rate : " + m_rate);
		System.out.println("The number of generation : " +generation);
		
		for (int i = 0; i <  generation; i++) {
			Population subPop = new Population(); // store new individuals are generated
			double[] popFitness = population.populationFitness();
			
			
			for (int j = 0; j < Population.populationLength / 2; j++) {
				Individual father = new  Individual();
				Individual mother = new  Individual();
				
				//use tournament to  select Individual 
				father = population.getIndividual(selection.touranmentSelection(popFitness, 2, MyRandom.r));
				mother = population.getIndividual(selection.touranmentSelection(popFitness, 2, MyRandom.r));		
				
					
				Individual offspring1 =  new Individual();
				Individual offspring2 =  new Individual();
			    
				// new  random number  [0;1]
				double d = 0 + (1 - 0) * MyRandom.r.nextDouble();
				// if random number >  crossover rate,  add  the individual to population without crossover 
		     	if (d > p_rate) {
					subPop.addIndiv(father);
					subPop.addIndiv(mother);
			     } else {
			//	 else do crossover, generate two offsprings, then  do mutation
			    	offspring1.setGene(crossover.clusterCrossover(father.getGene(), mother.getGene(),
	   		    		       num_vertex, ReadFiles.clusters, MyRandom.r));
			    	offspring2.setGene(crossover.clusterCrossover(father.getGene(), mother.getGene(),
		    		       num_vertex, ReadFiles.clusters,MyRandom.r)); 
//					offspring1.setGene(crossover.ClusterBFSCrossover(father.getGene(), mother.getGene(),
//   		    		       num_vertex, ReadFiles.clusters, MyRandom.r));
//					
//					offspring2.setGene(crossover.ClusterBFSCrossover(father.getGene(), mother.getGene(),
//	   		    		       num_vertex, ReadFiles.clusters, MyRandom.r));
					
					if( d < m_rate){
						 offspring1.setGene(mutation.mutationClusterTree(mother.getGene(),
								 num_vertex, ReadFiles.clusters, MyRandom.r));
						 offspring2.setGene(mutation.mutationClusterTree(father.getGene(), 
								 num_vertex, ReadFiles.clusters, MyRandom.r));
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
			Paint  p = new Paint();
			p.weightMatrix = population.getIndividual(bestFinessIndex).getGene();
			p.fitness = popFitness[bestFinessIndex];
//			p.fitness1 = popFitness1[bestFinessIndex];
			p.num_vertex = num_vertex;
		    gf.add(p);
		    gf.setVisible(true);
		    
		    
			// add the best Individual  of old Population to new Population
			subPop.population.set(MyRandom.r.nextInt(Population.populationLength), 
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
		double[][] temptable = 	population.getIndividual(bestFinessIndex).getGene();
		for(int tr = 0;tr < num_vertex; tr++ ){
			for(int p =0; p < num_vertex; p++){
				
				if(temptable[tr][p] !=temptable[p][tr]){
					System.out.println(tr+" "+p);
				}
				
				
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
//		 pw.println(popFitness[bestFinessIndex]);
		 
	}
//		  pw.close();
		  }
}

