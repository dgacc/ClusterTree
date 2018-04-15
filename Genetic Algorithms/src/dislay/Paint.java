package dislay;
import java.awt.*;
import java.util.ArrayList;
import java.util.Random;

import javax.swing.*;

import filesinout.ReadFiles;
import structures.Cluster;
import structures.Vertex;

public class Paint extends JPanel {
    public int num_vertex;
    private Vertex[] vertexs; 
    private  double[][] weightMatrix ;
    private ArrayList<Cluster> clusters;
    private double fitness;
    private int root;
    private int  x = 0;
    

    
    private  void doDrawing(Graphics g) {
    	int size = 10;
    	int numberOfCluster  = clusters.size();
        Graphics2D g2d = (Graphics2D) g;
        for( int i = 0; i < numberOfCluster; i ++ ){
        	int numberClusterVertex = clusters.get(i).getCluster().size();
        	Random r = new Random();
        	g2d.setColor(new Color(r.nextInt(256),r.nextInt(256),r.nextInt(256)));
        	for( int  j = 0; j < numberClusterVertex; j++ ){
        	g2d.fillOval((int)vertexs[clusters.get(i).getCluster().get(j)].getX()*size - 5 + x,
        			(int)vertexs[clusters.get(i).getCluster().get(j)].getY()*size - 5, 10, 10);
        	}
        }
        for(int i = 0; i < num_vertex; i++){
        	for( int j = 0; j < num_vertex; j++){
        		 g2d.setColor(Color.RED);
        		 g2d.drawString(""+( i), (int)vertexs[i].getX()*size + 10 + x, (int)vertexs[i].getY()*size + 10);
        		 if(weightMatrix[i][j] > 0){
        		 g2d.setColor(Color.BLUE);	 
                 g2d.drawLine((int)vertexs[i].getX()*size + x, (int)vertexs[i].getY()*size, 
                		 (int)vertexs[j].getX()*size + x, (int)vertexs[j].getY()*size);        		
        		 }
        }
   		}
        
        
        g2d.setColor(Color.BLACK);
        g2d.fillOval((int)vertexs[root].getX()*size - 10 + x ,
        		(int)vertexs[root].getY()*size - 10, 20, 20);
        g2d.drawString("fitness:" + fitness , 50 + x + 200, 20);
      
    }
 
    @Override
    public void paintComponent(Graphics g) {
 
        super.paintComponent(g);
        doDrawing(g);
    }
    
     public void  setPaint(double[][] weightMatrix , Vertex[] vertexs,ArrayList<Cluster> clusters,int num_vertex, 
    		 double fitness, int index, int root){
    	this.weightMatrix = weightMatrix;
    	this.vertexs = vertexs;
    	this.clusters = clusters;
    	this.num_vertex = num_vertex;
    	this.fitness = fitness;
    	this.x = index;
    	this.root = root;
    }
    
    
    
    
    
    
}
//
 
