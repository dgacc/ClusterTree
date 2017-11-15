package Cluster_Tree;
import java.awt.*;
import java.util.ArrayList;
import java.util.Random;

import javax.swing.*;

class Paint extends JPanel {
    public int num_vertex;
    private static Vertex[] vertexs = ReadFiles.vertices;
//    public double[][] weightMatrix ;
    public double[][] weightMatrix ;
    private static ArrayList<Cluster> clusters = ReadFiles.clusters;
    public double fitness;
    public double fitness1;
    
    
    
    
    private  void doDrawing(Graphics g) {
    	int x = 10;
    	int numberOfCluster  = clusters.size();
        Graphics2D g2d = (Graphics2D) g;
        for( int i = 0; i < numberOfCluster; i ++ ){
        	int numberClusterVertex = clusters.get(i).getCluster().size();
        	Random r = new Random();
        	g2d.setColor(new Color(r.nextInt(256),r.nextInt(256),r.nextInt(256)));
        	for( int  j = 0; j < numberClusterVertex; j++ ){
        	g2d.fillOval((int)vertexs[clusters.get(i).getCluster().get(j)].getX()*x - 5,
        			(int)vertexs[clusters.get(i).getCluster().get(j)].getY()*x - 5, 10, 10);
        	}
        }
        for(int i = 0; i < num_vertex; i++){
        	for( int j = 0; j < num_vertex; j++){
        		 g2d.setColor(Color.RED);
        		 g2d.drawString(""+( i), (int)vertexs[i].getX()*x + 10, (int)vertexs[i].getY()*x + 10);
        		 if(weightMatrix[i][j] > 0){
        		 g2d.setColor(Color.BLUE);	 
                 g2d.drawLine((int)vertexs[i].getX()*x, (int)vertexs[i].getY()*x, 
                		 (int)vertexs[j].getX()*x, (int)vertexs[j].getY()*x);        		
        		 }
        }
   		}
        
        
        g2d.setColor(Color.BLACK);
        g2d.fillOval((int)vertexs[ReadFiles.root].getX()*x - 10 ,
        		(int)vertexs[ReadFiles.root].getY()*x - 10, 20, 20);
        g2d.drawString("fitness:" + fitness , 50, 20);
        g2d.drawString("fitness1:" + fitness1 , 300, 20);
    }
 
    @Override
    public void paintComponent(Graphics g) {
 
        super.paintComponent(g);
        doDrawing(g);
    }
}
//
 
