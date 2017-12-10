/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
//import static javax.swing.Spring.scale;
//import static jogamp.nativewindow.SurfaceScaleUtils.scale;
import util.TFChangeListener;
import util.VectorMath;

import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private GradientVolume gradients = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    
    
    private boolean shading = false;  // color of the light
    private TFColor ambient = new TFColor();
    private TFColor diffuse = new TFColor();
    private int type = 1;
    private double k_ambient = 0.1;
    private double k_diffuse = 0.7;
    private double k_spec = 0.2;
    private double alpha1 = 20;
    private boolean Trilin = false;
 
    public void setShading()
    {
     if(shading == true){shading = false; this.changed();}
     else{
     shading = true; 
     this.changed();
     }
    }
    public void setTrilinear()
    {
     if(Trilin == true){Trilin = false; this.changed();}
     else{
     Trilin = true; 
     this.changed();
     }
    }

    // Return Phong shading model
    private TFColor getPhongShadedColor(double[] coord, TFColor voxelColor, double[] viewVec) {
        if (coord[0] < 0 || coord[0] >= volume.getDimX()
                || coord[1] < 0 || coord[1] >= volume.getDimY()
                || coord[2] < 0 || coord[2] >= volume.getDimZ()) {
            voxelColor.set(0, 0, 0, voxelColor.a);
            return voxelColor;
        }
        //Setting diffuse to voxelColor
        diffuse = voxelColor;
        VoxelGradient vg = gradients.getGradient(
                (int) Math.floor(coord[0]),
                (int) Math.floor(coord[1]),
                (int) Math.floor(coord[2]));
        //return getColor(viewVec, vg.getNormal(), viewVec, diffuse);
        //double[] L, double[] N, double[] V, TFColor dif
                
        double[] N = vg.getNormal();
        //reflection vector computation as given in the slides R = 2(N.L)N - L
        //Here L is viewVec
        double[] R = VectorMath.subtract((VectorMath.scale(N, 2*VectorMath.dotproduct(N, viewVec))), viewVec);
      
        //computing dot product L.N & (V.R)^alpha1
        //V is also viewVec
        double LN = Math.max(0, VectorMath.dotproduct(viewVec, N));
        double VR = VectorMath.dotproduct(viewVec, R);
        VR = (VR <= 0 || LN <= 0) ? 0 : Math.pow(VR, alpha1); 
        
        //Setting diffuse to white color
        ambient = new TFColor(1,1,1,1);
        //computing color using phong shading model
        TFColor shadedcolor = new TFColor(
                k_ambient * ambient.r + k_diffuse * diffuse.r * LN + k_spec * diffuse.r * VR,
                k_ambient * ambient.g + k_diffuse * diffuse.g * LN + k_spec * diffuse.g * VR,
                k_ambient * ambient.b + k_diffuse * diffuse.b * LN + k_spec * diffuse.b * VR,
                diffuse.a);
        return shadedcolor;
    }
    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;

        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        // create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases
        // linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        
        // uncomment this to initialize the TF with good starting values for the orange dataset 
        tFunc.setTestFunc();
        
        
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }
    double getVoxel(double[] coord) {
        if(!Trilin){
                if (coord[0] < 0 || coord[0] > volume.getDimX() || coord[1] < 0 || coord[1] > volume.getDimY()
                || coord[2] < 0 || coord[2] > volume.getDimZ()) {
            return 0;
        }

        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);

        return volume.getVoxel(x, y, z);
        }
        else{
        
            double x = coord[0];
            double y = coord[1];
            double z = coord[2];

            int x0 = (int) Math.floor(x);
            int x1 = (int) Math.min(Math.ceil(x), volume.getDimX() - 1);
            int y0 = (int) Math.floor(y);
            int y1 = (int) Math.min(Math.ceil(y), volume.getDimY() - 1);
            int z0 = (int) Math.floor(z);
            int z1 = (int) Math.min(Math.ceil(z), volume.getDimZ() - 1);

             if (coord[0] < 0 || coord[0] > volume.getDimX() || coord[1] < 0 || coord[1] > volume.getDimY()
                    || coord[2] < 0 || coord[2] > volume.getDimZ() || x0 >= volume.getDimX() || x0 < 0 || y0 >= volume.getDimY() || y0 < 0 || z0 >= volume.getDimZ() || z0<0 || x1<0 || x1 > volume.getDimX() || y1 < 0 || y1 > volume.getDimY() || z1<0 || z1> volume.getDimZ()) {
                return 0;
            }

            // Getting the values of the 8 points
            double v000 = volume.getVoxel(x0, y0, z0);
            double v001 = volume.getVoxel(x0, y0, z1);
            double v010 = volume.getVoxel(x0, y1, z0);
            double v100 = volume.getVoxel(x1, y0, z0);
            double v101 = volume.getVoxel(x1, y0, z1);
            double v011 = volume.getVoxel(x0, y1, z1);
            double v110 = volume.getVoxel(x1, y1, z0);
            double v111 = volume.getVoxel(x1, y1, z1);

            // interpolating about x dimension
            double v00 = ((1-((x-x0)/(x1-x0)))*v000 )+ (((x-x0)/(x1-x0))*v100);
            double v01 = ((1-((x-x0)/(x1-x0)))*v001 )+ (((x-x0)/(x1-x0))*v101);
            double v10 = ((1-((x-x0)/(x1-x0)))*v010 )+ (((x-x0)/(x1-x0))*v110);
            double v11 = ((1-((x-x0)/(x1-x0)))*v011 )+ (((x-x0)/(x1-x0))*v111);

            // interpolating about y dimension
            double v0 = ((1-((y-y0)/(y1-y0)))*v00 )+ (((y-y0)/(y1-y0))*v10);
            double v1 = ((1-((y-y0)/(y1-y0)))*v01 )+ (((y-y0)/(y1-y0))*v11);

            // interpolating about z dimension
            double v = ((1-((z-z0)/(z1-z0)))*v0 )+ (((z-z0)/(z1-z0))*v1);
            return v;
        }
    }  

    void slicer(double[] viewMatrix) {
        int res_param = 1;
        if(interactiveMode) 
        {
            res_param = 5;
        }
        
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        
        for (int j = 0; j < image.getHeight(); j=j+res_param) {
            for (int i = 0; i < image.getWidth(); i=i+res_param) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];

                double val = getVoxel(pixelCoord);
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                //voxelColor = tFunc.getColor(val);
                
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                for (int m = 0; m <res_param; m++) {
                    for (int n= 0; n <res_param; n++) {
                        if (n + i < image.getHeight() 
                            && n + i >= 0
                            && m + j < image.getWidth() 
                            && m + j >= 0) {
                            image.setRGB(n + i, m + j,pixelColor);
                        }
                    }
                }
            }
        }

    }
    
      
    void mip(double[] viewMatrix) {
        int res_param = 1;
        if (interactiveMode) {
            res_param = 3  ;
        }
        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }
        

        double[] uVec = new double[3];
        double[] vVec = new double[3];
        double[] viewVec = new double[3];     
        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        //image is square
        int imageCenter = image.getWidth() / 2;
        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        //sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();
        for (int j = 0; j < image.getHeight(); j+= res_param) {
            for (int i = 0; i < image.getWidth(); i+= res_param) {
                
                int diag = (int) Math.sqrt(volume.getDimX() * volume.getDimX() + volume.getDimY() * volume.getDimY() + volume.getDimZ() * volume.getDimZ()); 
                double maxIntensity = 0;
                 
                //taking along the ray axis , from 0 to diagonal length we have taken as diag gives max length and we dont miss the data 
                for (int step_k = 0; step_k < diag - 1; step_k++) {

                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) +(step_k - imageCenter) * viewVec[0] + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + (step_k - imageCenter) * viewVec[1] + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) +  (step_k - imageCenter)*viewVec[2] + volumeCenter[2];


                    double val = getVoxel(pixelCoord);
                
                
                    if (val > maxIntensity) {
                        maxIntensity = val;
                    }
                    if (val / max > 0.95) {
                        break;
                    }
                }
               
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = maxIntensity/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = maxIntensity > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                //voxelColor = tFunc.getColor(val);
                
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                
                for (int m = 0; m <res_param; m++) {
                    for (int n= 0; n <res_param; n++) {
                        if (n + i < image.getHeight() 
                            && n + i >= 0
                            && m + j < image.getWidth() 
                            && m + j >= 0) {
                            image.setRGB(n + i, m + j,pixelColor);
                        }
                    }
                }
            }
        }

    }
    private void composite(double[] viewMatrix) {
        int res_param = 1;
        if(interactiveMode) {
            res_param = 5;
        }
        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        double[] viewVec = new double[3];
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.normalize(viewVec);
        

        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        TFColor voxelColor, voxelColor1 = new TFColor();
    
        for (int j = 0; j < image.getHeight(); j=j+res_param) {
            for (int i = 0; i < image.getWidth(); i=i+res_param) {

                voxelColor1.set(0, 0, 0, 1);
                //double[] v = VectorMath.subtract(farP, nearP);
                int diag = (int) Math.sqrt(volume.getDimX() * volume.getDimX() + volume.getDimY() * volume.getDimY() + volume.getDimZ() * volume.getDimZ()); 
                
                for (int step_k =0;step_k < diag - 1;step_k = step_k + 1) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) 
                        + vVec[0] * (j - imageCenter) 
                         +  viewVec[0]*(step_k - imageCenter ) + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter)
                        + vVec[1] * (j - imageCenter) 
                         + viewVec[1]*(step_k - imageCenter) + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) 
                        + vVec[2] * (j - imageCenter)
                         + viewVec[2]*(step_k - imageCenter)+volumeCenter[2];

                double val = getVoxel(pixelCoord);

                // apply the transfer function to obtain a color
                voxelColor = tFunc.getColor(val);

                //phong shading
                    if (shading) {
                        voxelColor = getPhongShadedColor(pixelCoord, voxelColor, viewVec);
                    }
                    
                    //compositing
                    voxelColor1.r = voxelColor.a * voxelColor.r + (1 - voxelColor.a) * voxelColor1.r;
                    voxelColor1.g = voxelColor.a * voxelColor.g + (1 - voxelColor.a) * voxelColor1.g;
                    voxelColor1.b = voxelColor.a * voxelColor.b + (1 - voxelColor.a) * voxelColor1.b;
                }
 
                //BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor1.a <= 1.0 ? (int) Math.floor(voxelColor1.a * 255) : 255;
                int c_red = voxelColor1.r <= 1.0 ? (int) Math.floor(voxelColor1.r * 255) : 255;
                int c_green = voxelColor1.g <= 1.0 ? (int) Math.floor(voxelColor1.g * 255) : 255;
                int c_blue = voxelColor1.b <= 1.0 ? (int) Math.floor(voxelColor1.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                for (int m = 0; m <res_param; m++) {
                    for (int n= 0; n <res_param; n++) {
                        if (n + i < image.getHeight() 
                            && n + i >= 0
                            && m + j < image.getWidth() 
                            && m + j >= 0) {
                            image.setRGB(n + i, m + j,pixelColor);
                        }
                    }
                }
            }
        }
    }
     void Transfer2DFunction (double[] viewMatrix){
         
        int res_param = 1;
        if (interactiveMode) {
            res_param = 5 ;
        }
          // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        
                // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];  
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);

        // image is square
        int imageCenter = image.getWidth() / 2;
        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();

        // Get the maximum length that is diagonal
        int diag = (int) Math.sqrt(volume.getDimX() * volume.getDimX() + volume.getDimY() * volume.getDimY() + volume.getDimZ() * volume.getDimZ());

        // For each pixel of the image
        for (int j = 0; j < image.getHeight(); j+=res_param) {
            for (int i = 0; i < image.getWidth(); i+=res_param) {
                // First, set a color variable in which we can put all the colors together
                TFColor voxelColor = new TFColor(0,0,0,1);

                for(int step_k = 0; step_k < diag - 1; step_k++) {
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter) + viewVec[0] * (step_k - imageCenter) + volumeCenter[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter) + viewVec[1] * (step_k - imageCenter) + volumeCenter[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter) + viewVec[2] * (step_k - imageCenter) + volumeCenter[2];
                  
                    int val = (int) getVoxel(pixelCoord);
                    if (val == 0) {
                        continue;
                    }
                    TransferFunction2DEditor.TriangleWidget triagWid = this.tfEditor2D.triangleWidget;
                    float baseIntensity = triagWid.baseIntensity;
                    TFColor color = triagWid.color;
                    double radius = triagWid.radius;
                    double maxGradient = triagWid.upGradient;
                    double minGradient = triagWid.downGradient;
                    double opacity = 0;
                    
                    VoxelGradient gradient = gradients.getGradient((int) pixelCoord[0], (int) pixelCoord[1], (int) pixelCoord[2]);
                    if (gradient.mag > maxGradient || gradient.mag < minGradient){
                        opacity=0;
                    }
                    else if (val == baseIntensity && gradient.mag == 0){
                        opacity = 1;
                    } else if (gradient.mag > 0 && ((val - (radius * gradient.mag)) <= baseIntensity && baseIntensity <= (val + (radius * gradient.mag)))) {
                            opacity = color.a * (1 - (1 / radius) * Math.abs(((baseIntensity - val) / gradient.mag)));
                        } else {
                        opacity = 0;
                    }
                    
                    
                    if (shading) {
                        color = getPhongShadedColor(pixelCoord, color,viewVec);
                    } 
                    
                    voxelColor.r = (color.r * opacity ) + (voxelColor.r * (1 - opacity));
                    voxelColor.g = (color.g * opacity ) + (voxelColor.g * (1 - opacity));
                    voxelColor.b = (color.b * opacity ) + (voxelColor.b * (1 - opacity)); 
                  }
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
                
               for (int m = 0; m <res_param; m++) {
                    for (int n= 0; n <res_param; n++) {
                        if (n + i < image.getHeight() 
                            && n + i >= 0
                            && m + j < image.getWidth() 
                            && m + j >= 0) {
                            image.setRGB(n + i, m + j,pixelColor);
                        }
                    }
                }
            }
        }
     }
  
    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }

    @Override
    public void visualize(GL2 gl) {


        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
            
            switch (type){
            case 1:
                slicer(viewMatrix);
                break;
            case 2:
                mip(viewMatrix);
                break;
            case 3:
                composite(viewMatrix);
                break;
            case 4:
                Transfer2DFunction(viewMatrix);
                break;
            default:
                slicer(viewMatrix);
        }    
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];
    public void setType(int type) {
        this.type = type;
        this.changed();
    }
    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}
