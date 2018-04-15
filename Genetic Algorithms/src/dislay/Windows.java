package dislay;

import javax.swing.*;

public class Windows  extends JFrame{
	JLayeredPane pane = getLayeredPane();
	
    private JFrame gf = new JFrame();
    public void runWindow(String title){
		gf.setVisible(true);
		gf.setSize(1400, 800);
		gf.setTitle(title);
		gf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		gf.setVisible(true);
//		pane.add(gf);
//		gf.setLayout(null);

	}	
    
    
     public  void addPaint(Paint paint){
    	 gf.add(paint);
    	 gf.setVisible(true);

     }
    
}
