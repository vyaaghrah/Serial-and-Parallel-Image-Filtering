import java.io.File;
import java.io.IOException;
import java.awt.image.BufferedImage;
import java.util.Arrays;
import javax.imageio.ImageIO;
public class Medianfilter_Serial 
{
    public static void main(String args[])
    {
        BufferedImage img=null;
        File f=null;
        
        long startTime = System.currentTimeMillis();
        try{
            f=new File("/home/addy/Desktop/HPC/Coins.jpeg");
            img=ImageIO.read(f);
        }
	catch(IOException e)
        {
            System.out.println(e);
        }
        
        int width=img.getWidth();
        int height=img.getHeight();
        
        for(int y=1;y<height-1;y++)
            for(int x=1;x<width-1;x++)
            {
                int ar[]=new int[9];
                ar[0]=img.getRGB(x, y);
                ar[1]=img.getRGB(x-1, y-1);
                ar[2]=img.getRGB(x-1, y);
                ar[3]=img.getRGB(x-1, y+1);
                ar[4]=img.getRGB(x, y-1);
                ar[5]=img.getRGB(x, y+1);
                ar[6]=img.getRGB(x+1, y-1);
                ar[7]=img.getRGB(x+1, y);
                ar[8]=img.getRGB(x+1, y+1);
                
		int i, j;
	        for (i = 0; i < 8; i++)       
		for (j = 0; j < 8-i; j++)  
        	if (ar[j] > ar[j+1]) 
              	{
			int temp=ar[j];
			ar[j]=ar[j+1];
			ar[j+1]=temp;
		}

                //Arrays.sort(ar);
                int p=img.getRGB(x, y);
                int p1=ar[4];
                
                if(p>255||p<0)
                    img.setRGB(x, y, p1);
            }
        
        try
	{
	      ImageIO.write(img, "bmp", new File("./CoinsSerial.bmp"));
        }
	catch(IOException e)
        {
            System.out.println(e);
        }
        long endTime = System.currentTimeMillis();
        System.out.println(Math.abs(startTime-endTime)+" milliseconds");
    }
}
