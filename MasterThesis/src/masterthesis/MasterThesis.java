/**
 * This package contains various methods used during my master thesis.
 * The main steps were the following: data filtration to reduce noisiness,
 * breaking up the raw data file into smaller segments (a year's worth of data
 * adds up to huge files that are hard to handle in Java and Excel), and model
 * fitting. The method of main interest is predict_exponentials(), where an
 * outside program (Matlab) is called and its output from the command line is
 * read back into this script.
 * 
 * As a one-time script collector, I must apologize for the structure;
 * I tried my best to make it more readable.
 * 
 * Author: Mátyás Tóth, code written: 2016-2017 in NetBeans IDE 8.0.2.
 */
package masterthesis;

import java.awt.Toolkit;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import utilsforgeneraluse.UtilsForGeneralUse;

public class MasterThesis {

    //Input file
    static File file;
    
    //Output file
    static File outfile;
    
    //Data field separator.
    static String DS = ";";
    
    static double bound;
    
    //The Y column for 2D plotting.
    static int ycolumn = -1;
    
    //The interval used for filtering - further description at the method.
    static int interval = -1;
    static double standard_deviation = 0;
    static double limiting_standard_deviation = 0.001;
    
    
    /*
     * The main method is used only as a setup to run the desired script.
     */
    public static void main(String[] args) throws IOException, Exception {
        //densviscosity();
        //separator();
        
        file = new File("file_location");
        outfile = new File("output_file_location");

        doAverage(0.1, 12, 4);
        file = new File("another_file_location");
        predict(file, 60.1425, false);
        //predict_exponentials(2.42, 24.8);
    }

    /**
     * This method was used to fit two different exponential functions to the
     * filtered data. In order to simulate a real-time experiment, it runs
     * through the data file recursively; first the models are fitted on the
     * first five data points, then on the first six, and so on.
     * 
     * In order to spare errors and effort, this script makes calls through the
     * command line to Matlab to fit the exponential models. After each fit,
     * the results are read back to this script and written to the specified
     * output file.
     * 
     * @param Rlimit The limiting membrane resistance which the exponentials
     *                  are infinitely getting close to.
     * @param time  The end time of the run, given by the plant.
     *              At this point, the membranes were chemically cleaned.
     * @throws IOException 
     */
    public static void predict_exponentials(double Rlimit, double time) throws IOException {
        /*Rows: different parameters, first column: fitted value,
         second: lower error bound, third: upper error bound
         */
        double[][] parameters = null;
        double prediction;

        //File regarded as input.
        file = new File("file:\\...\\input.csv");

        //Files containing results of fitting.
        File outfile_exp = new File("file:\\...\\matlabfit_exp.csv");
        File outfile_strexp = new File("file:\\...\\matlabfit_strexp.csv");

        //Temporary file, containing limited amount of data;
        //  it's regarded as 'this much is known at this point'.
        File temp = new File("file:\\...\\temp.csv");
        if (temp.exists()) {
            temp.delete();
        }

        //Script location for exponential fitting
        File script1 = new File("file:\\...\\exp_fit.m");

        //Script location for stretched exponential fitting
        File script2 = new File("file:\\...\\strexp_fit.m");

        BufferedReader br = new BufferedReader(new FileReader(file));
        BufferedWriter bw_exp = new BufferedWriter(new FileWriter(outfile_exp));
        BufferedWriter bw_strexp = new BufferedWriter(new FileWriter(outfile_strexp));
        BufferedWriter bw;

        //Fitting needs a few initial points; I used five in this example.
        int maxlines = 5;
        
        //This shows how many runs have we completed.
        int pointer;
        
        //There, I am making a reference to an outside project; it's only
        // used to determine how many lines are present in a given file.
        int filelines = UtilsForGeneralUse.countLine(file);
        
        //An array containing the time stamps...
        double[] t = new double[filelines];
        for (int i = 0; i < t.length; i++) {
            //...read from the first column of an input file.
            t[i] = Double.valueOf(br.readLine().split(DS)[0]);
        }
        br.close();

        String line = "";
        //If the time at which the limiting resistance is approached closely,
        //  exit the procedure.
        boolean EXIT = false;
        
        //Iterate through the file.
        while (maxlines < filelines + 1) {
            //If by any chance, the temporary file still exists, remove it.
            if (!temp.exists()) {
                while (true) {
                    try {
                        temp.createNewFile();
                        break;
                    } catch (IOException ioe) {
                        //In case of an error, give a sign to the unaware!
                        Toolkit.getDefaultToolkit().beep();
                    }
                }
            }
            br = new BufferedReader(new FileReader(file));
            bw = new BufferedWriter(new FileWriter(temp));

            //Writing temporary input file for fitting.
            pointer = 0;
            while (pointer < maxlines) {
                line = br.readLine();
                //This makes sure that irrevelant values are not included in
                //  the fitting runs.
                if (time < Double.valueOf(line.split(DS)[0])) {
                    EXIT = true;
                }
                //Another method from an outside project.
                //Prints a line with BufferedWriter.
                UtilsForGeneralUse.println(bw, line);
                //Making sure that the next fitting run goes one data point further.
                pointer++;
            }
            br.close();
            bw.close();

            if (EXIT) {
                break;
            }

            //Print which time point are we at, and what is the end point.
            //This is mostly a feedback to know how the script is progressing.
            System.out.println(Double.valueOf(line.split(DS)[0]) + "\t" + time);

            //Exponential fit to given file, and getting back the fitting parameters.
            parameters = matlab(script1, 0);
            //Get the expected time at which the membrane is fouled.
            prediction = (-1) * Math.log(1 - (Rlimit / parameters[0][0])) / parameters[1][0];
            //Print everything else too.
            UtilsForGeneralUse.print(bw_exp, Double.toString(t[maxlines - 1]) + DS);
            UtilsForGeneralUse.print(bw_exp, Double.toString(prediction) + DS);
            for (int i = 0; i < parameters.length; i++) {
                for (int j = 0; j < parameters[0].length; j++) {
                    if (parameters[i][j] != -1) {
                        UtilsForGeneralUse.print(bw_exp, Double.toString(parameters[i][j]) + DS);
                    }
                }
            }
            UtilsForGeneralUse.println(bw_exp, "");
            
            //The same procedure, but with the stretched exponential model.
            //Stretched exponential fit to given file.
            parameters = matlab(script2, 1);
            prediction = Math.pow((-1) * Math.log(1 - (Rlimit / parameters[0][0])), 1 / parameters[2][0]) / parameters[1][0];
            //Print results...
            UtilsForGeneralUse.print(bw_strexp, Double.toString(t[maxlines - 1]) + DS);
            UtilsForGeneralUse.print(bw_strexp, Double.toString(prediction) + DS);
            for (int i = 0; i < parameters.length; i++) {
                for (int j = 0; j < parameters[0].length; j++) {
                    if (parameters[i][j] != -1) {
                        UtilsForGeneralUse.print(bw_strexp, Double.toString(parameters[i][j]) + DS);
                    }
                }
            }
            UtilsForGeneralUse.println(bw_strexp, "");

            temp.delete();
            maxlines++;
        }

        if (temp.exists()) {
            temp.delete();
        }

        bw_exp.close();
        bw_strexp.close();
    }

    /**
     *
     * @param input File to be used as input.
     * @param mode '0' for exponential fit, '1' for stretched exponential fit.
     * @return Double array of fitted parameters. Null if error occurred.
     *          Optionally, a double array filled with -1 if the error is not
     *          crashing the matlab script.
     * @throws IOException
     */
    public static double[][] matlab(File input, int mode) throws IOException {
        //Single exponential.
        if (mode != 0 || mode != 1) {
            System.err.println("Mode not supported!");
            return null;
        }
        double[][] rd = new double[mode + 2 + 5][3];
        for (int i = 0; i < rd.length; i++) {
            //Initialize the array and fill it with -1.
            Arrays.fill(rd[i] = new double[rd[0].length], -1);
        }

        //Initializing a runtime process which the JVM can follow.
        Runtime rt = Runtime.getRuntime();
        
        //Starting the process with a special script to keep Matlab minimized
        //  and force it to output the results to the command line.
        
        //As a note, this mode of Matlab was intended to be used as a
        //  debugging mode - thus the output is an error stream.
        Process proc = rt.exec("matlab -nosplash -nodesktop -nojvm -nodisplay -minimize -r \"run('" + input.getAbsolutePath() + "')\"  -wait -log");
        BufferedReader stdError = new BufferedReader(new InputStreamReader(proc.getErrorStream()));
        
        
        String s = null;
        //Reading the output of the process.
        while ((s = stdError.readLine()) != null) {
            
            //When we've got to the coefficients of the fitting...
            if (s.contains("Coefficients")) {
                int i = 0;
                int j = 0;
                //As these coefficients can be separated by empty lines, this
                //  cycle is necessary.
                while ((s = stdError.readLine()) != null && !s.isEmpty()) {
                    //Reformatting the line read; removing unnecessary
                    //  characters and turning the decimal '.' into ','.
                    String sh = s.replace(" ", "").replace("Rss", "").replace("a", "")
                            .replace("b", "").replace("=", "").replace("(", ",").replace(")", "");
                    j = 0;
                    //Now all numbers are separated by comma.
                    // Reading them into the output array...
                    for (String t : sh.split(",")) {
                        try {
                            rd[i][j] = Double.parseDouble(t);
                        } catch (Exception e) {
                            return null;
                        }
                        j++;
                    }
                    i++;
                }
                
                /*
                    Apologies for my language.
                    Anyway, fetching the std. error of the fitting.
                */
                while ((s = stdError.readLine()) != null && !s.contentEquals("done mothafucka")) {
                    if (s.isEmpty() || s.contains("gof")) {
                        continue;
                    }
                    String sh = s.substring(16);
                    //System.out.println(sh);
                    rd[i][0] = Double.parseDouble(sh);
                    i++;
                }
            }
            
            //If the fitting was unsuccessful, fill the array with '-1', and
            //  return. This included a bit more precise error message and
            //  tracking, but after a few rounds, it was no longer revelant.
            if (s.contains("Error")) {
                System.err.println("OH NOES");
                Arrays.fill(rd, -1);
                return rd;
            }
        }
        return rd;
    }

    /**Deprecated method. It did a simple averaging on columns.
     * 
     * @deprecated 
     * @param diff Difference of maximum and minimum value between which
     * averaging takes place.
     * @param column_X
     * @param column_Y
     * @throws IOException
     */
    public static void doAverage(double diff, int column_X, int column_Y) throws IOException {
        BufferedReader br = new BufferedReader(new FileReader(file));
        BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));

        String line;
        String[] tokens;
        int linecount = 0;

        while ((line = br.readLine()) != null) {
            linecount++;
        }
        br.close();

        br = new BufferedReader(new FileReader(file));
        double[][] array = new double[linecount][2];
        linecount = 0;
        while ((line = br.readLine()) != null) {
            tokens = line.split(DS);
            array[linecount][0] = Double.parseDouble(tokens[column_X]);
            array[linecount][1] = Double.parseDouble(tokens[column_Y]);
            linecount++;
        }
        br.close();

        double[] av = new double[linecount];

        int j = 0;
        double sum = 0;
        for (int i = 0; i < av.length; i++) {
            j = i;
            sum = 0;
            while ((diff > array[i][0] - array[j][0]) && j > 0) {
                sum += array[j][1];
                j--;
            }
            av[i] = sum / (i - j);
        }

        br = new BufferedReader(new FileReader(file));
        int i = 0;
        while ((line = br.readLine()) != null) {
            bw.write(line + DS + av[i]);
            i++;
            bw.newLine();
            bw.flush();
        }
        br.close();
        bw.close();
    }

    /**
     * This method filters out data points using given conditions. The actual
     * filtering conditions are given in
     * {@link masterthesis.MasterThesis#isWithinBounds(double[], double, int, int)}.
     *
     * @param f Input file.
     * @param outf Output file.
     * @param interval Number of values used for calculation of averages. Should
     * be an uneven integer, at least 1. (Preferably 5 or 7.)
     * @param limit_sd shows the maximum standard deviation of the data
     * accepted. This is chosen arbitrarily - there was a lot of data to cull!
     * @param ycolumn is the number of column used (if there is more than one)
     * for the filtering method.
     * @param isTest Tells the method whether this is a test run or not.
     * @return A @double value is returned, which is the % of the passed vs. all
     * data.
     * @throws IOException In case something weird happens, like no file can be
     * written/read.
     */
    public static double filter(File f, File outf, int interval, double limit_sd, int ycolumn, boolean isTest) throws IOException {
        BufferedWriter bw = null;
        
        int percentage = 0;
        int count = 0;
        
        //Try reading the file...
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            //To be honest, I do not know why I kept a testing case in the code.
            if (!isTest) {
                bw = new BufferedWriter(new FileWriter(outfile));
            }
            //Initializing the line and the field tokens...
            String line = "";
            String[] tokens;
            
            //The number and string sequence containing the numbers from a line.
            // Strictly speaking, this could have been done better, and 'lineseq'
            // is not necessary.
            double[] sequence = new double[interval * 2 + 1];
            String[] lineseq = new String[interval * 2 + 1];
            int start = 0;
            //The index of the number that needs to be investigated - deciding 
            //  whether it belongs to the data set or not.
            int ind = (start + Math.floorDiv(sequence.length, 2)) % sequence.length;
            //Start reading the values from the original log file.
            //  NOTE: this is only reading the first few lines to start the
            //        algorithm.
            for (int i = 0; i < sequence.length; i++) {
                line = br.readLine();
                lineseq[i] = line;
                sequence[i] = Double.valueOf(line.split(DS)[ycolumn]);
            }
            //If the inspected value is within bounds, write it to the output file.
            if (isWithinBounds(sequence, ind, interval, limit_sd)) {
                if (!isTest) {
                    bw.write(lineseq[ind]);
                    bw.newLine();
                    bw.flush();
                }
                //Statistical stuff; show how much data was culled.
                percentage++;
            }
            count++;
            
            //Continuing with the rest of the file, now that the array containing
            //  the first few values is initialized.
            while ((line = br.readLine()) != null) {
                /*
                  This is a neat trick to avoid push() and pop() from ArrayLists,
                  but a bit harder to explain.
                  Also, this is unnecessarily complicating it compared to ArrayLists.
                  
                  This cycle overwrites the 'oldest' data in the array, and
                  does the evaluation on the data following the previous
                  investigated value. As such, it 'loops' through the array, as
                  if it is shifted to the left by one, and a number is added to
                  the right.
                */
                
                //Otherwise, does perfectly the same as before.
                start++;
                ind++;
                start %= sequence.length;
                ind %= sequence.length;

                lineseq[start] = line;
                tokens = line.split(DS);
                sequence[start] = Double.valueOf(tokens[ycolumn]);

                if (isWithinBounds(sequence, ind, interval, limit_sd)) {
                    if (!isTest) {
                        bw.write(lineseq[ind]);
                        bw.newLine();
                        bw.flush();
                    }
                    percentage++;
                }
                count++;
            }
        } catch (Exception e) {
            //A bit more detailed error message.
            System.err.println(e.getCause());
            return -1;
        }
        if (bw != null) {
            bw.close();
        }
        
        //Some statistics; reported in project, but has little meaning.
        return (100 * (double) percentage / (double) count);
    }

    /**
     * This method checks whether the value in the middle of the current array
     * is within error bounds of the average of previous/neighbouring/following
     * values or not.
     * 
     * The array is presumed to contain an uneven number of doubles. In case of
     * an array containing 21 numbers the procedure looks as follows:
     * 
     * n1 n2 n3 ... V ... nx, ny, nz
     *   where V is the value to be inspected.
     * 
     * (Note that V can be anywhere in the array, not necessarily in the middle.)
     * 
     * This method runs a series of checks, where first the std. deviance of the
     * first N numbers (n1, n2, n3...) is calculated and compared to V.
     * V is always skipped.
     * Then the same procedure is done, but starting from n2, n3, ... until the
     * last number is nz.
     * 
     * If V is within the maximum allowed deviance (limit_sd) at least once,
     * the program flags it as 'passed', true.
     *
     * @param array containing a (positive) number of variables, all of them
     * valid double values.
     * @param inv_ind The index of the investigated value.
     * @param N The amount of numbers taken into account when computing
     * averages.
     * @param limit_sd A double typically not higher 1; shows how close V has to
     *                  be to the standard deviation of the N numbers.
     * @return {@code true}, if at least one of the deviations is smaller than
     * the maximal allowed deviation.
     */
    public static boolean isWithinBounds(double[] array, int inv_ind, int N, double limit_sd) {
        
        //Average.
        double av = 0;
        double variance = 0;
        //Sum of deviances from the average.
        double sumdev = 0;
        double sd = 0;
        //Boolean showing whether there was a negative number in the N numbers.
        boolean nonegative = true;
        //Integer showing which n number should the script start from.
        int subindex;
        
        //Boolean indicating whether V passed.
        boolean wastrue = false;
        //Note that this is counting backwards!
        // We are iterating through N+1 values, and we simply check whether
        //  the indexed number is V.
        for (int index = N + 1; index > -1; index--) {
            //Start array from this point.
            subindex = (array.length + inv_ind - index) % array.length;

            for (int i = 0; i < N + 1; i++) {
                if (inv_ind != (array.length + subindex + i) % array.length) {
                    if (array[(array.length + subindex + i) % array.length] < 0) {
                        nonegative = false;
                        break;
                    }
                    av += array[(array.length + subindex + i) % array.length];
                }
            }
            av /= N;

            if (nonegative) {
                //Now standard deviation calculation.
                for (int i = 0; i < N + 1; i++) {
                    if (inv_ind != (array.length + subindex + i) % array.length) {
                        sumdev += Math.pow(av - array[(array.length + subindex + i) % array.length], 2);
                    }
                }
                variance = sumdev / N;
                sd = Math.sqrt(variance);
                
                //If the average of the values (negated by V) is smaller than
                //  the standard deviation set as a limit, then report it as true.
                if (sd * limit_sd > Math.abs(av - array[inv_ind])) {
                    wastrue = true;
                    standard_deviation = sd;
                    break;
                }
            }
            //Resetting values.
            av = 0;
            nonegative = true;
            sumdev = 0;
            variance = 0;
            sd = 0;
        }
        return wastrue;
    }

    /**
     * This method splits the raw data file into smaller segments, each one
     * starts from the backflush part of the cycle.
     * 
     * As such, there is not much interesting going on in this method.
     *
     * @throws IOException
     */
    public static void separator() throws IOException {
        file = new File("C:\\Users\\Matt\\university\\Aalborg\\master\\data\\paper_relaxations_data\\data.csv");
        BufferedReader br = new BufferedReader(new FileReader(file));
        File dir = new File("C:\\Users\\Matt\\university\\Aalborg\\master\\data\\paper_relaxations_data\\records");

        if (!dir.exists()) {
            dir.mkdir();
        }

        String line;
        ArrayList<String> records = new ArrayList<>();
        line = br.readLine();
        while (line != null) {

            if (Double.valueOf(line.split(DS)[1]) < 0) {
                records.add(line);
            }
            while ((line = br.readLine()) != null && Double.valueOf(line.split(DS)[1]) < 0) {
                records.add(line);
            }
            while ((line = br.readLine()) != null && Double.valueOf(line.split(DS)[1]) > 0) {
                records.add(line);
            }

            if (Double.valueOf(records.get(0).split(DS)[1]) < 0 && line != null) {
                if (outfile != null && outfile.getAbsolutePath().contentEquals("C:\\outfile" + records.get(0).split(DS)[0] + ".csv")) {
                    System.err.println(outfile.getAbsolutePath() + ":\tSameNameException!");
                    return;
                }
                outfile = new File("C:\\Users\\Matt\\university\\Aalborg\\master\\data\\paper_relaxations_data\\records\\"
                        + records.get(0).split(DS)[0] + ".csv");
                try (BufferedWriter bw = new BufferedWriter(new FileWriter(outfile))) {
                    for (String record : records) {
                        bw.write(record);
                        bw.newLine();
                        bw.flush();
                    }
                }
            }
            //Cleaning up, although this is unnecessary in Java.
            records.clear();
        }

    }

    /**
     * This method was used to extract and calculate some important physical variables
     * from the already available data.
     * 
     * While mostly simple, it contains a 5 point Lagrange interpolation to
     * estimate the density of the sewage water; mainly because the densities 
     * were measured with a minor delay.
     * @throws IOException 
     */
    public static void densviscosity() throws IOException {
        file = new File("C:\\input_file1");
        BufferedReader br = new BufferedReader(new FileReader(file));
        String line = "";

        String[] tokens;
        double[] dens = new double[41];
        while ((line = br.readLine()) != null) {
            tokens = line.split(DS);
            dens[Integer.valueOf(tokens[0])] = Double.valueOf(tokens[1]);
        }

        br.close();
        file = new File("C:\\input_file2");
        outfile = new File("C:\\output_file");
        br = new BufferedReader(new FileReader(file));
        BufferedWriter bw = new BufferedWriter(new FileWriter(outfile));

        DS = ";";
        double x;
        int j;
        //The following three variables are used for the Lagrange interpolation.
        double p;
        double p2m1;
        double p2m4;
        
        //...other variables necessary for the thesis.
        double density;
        double visc;
        double g = 9.81;
        //Trans-membrane pressure.
        double TMP;
        //Total resistance of the membrane.
        double Rtot;

        //Unit : micropascal. Base viscosity at 20 °C.
        double visc20 = 1002;
        while ((line = br.readLine()) != null) {
            tokens = line.split(DS, 8);
            x = Double.valueOf(tokens[4]);
            j = (int) Math.floor(x);
            p = x - j;
            p2m1 = p * p - 1;
            p2m4 = p2m1 - 3;

            //And 5 point Lagrange interpolation to get density!
            density = p * (p - 2) * p2m1 * dens[j - 2] / 24
                    - p * (p - 1) * p2m4 * dens[j - 1] / 6
                    + p2m1 * p2m4 * dens[j] / 4
                    - p * (p + 1) * p2m4 * dens[j + 1] / 6
                    + p * (p + 2) * p2m1 * dens[j + 2] / 24;

            //Now calculating viscosity; this is a written out mathematical formulae.
            visc = visc20 * Math.pow(10,
                    (20 - x) / (x + 96) * (1.2364 - 0.00137 * (20 - x) + 0.0000057 * (20 - x) * (20 - x)));

            TMP = g * density * Double.valueOf(tokens[7]) / 100;
            
            //Physical constants.
            int A = 40;
            int t = 3600;
            //Times 1 000 000, as the viscosity is given in MICROpascal.
            Rtot = TMP / (visc * Double.valueOf(tokens[1])) * 1000000 * A * t;
            
            //This is done to correct all SI units to the most used ones.
            //Now, Rtot has the "1/m" unit.

            double time = (Double.valueOf(tokens[0])) * 3600 * 24;
            bw.write(time + DS);
            for (int i = 1; i < tokens.length; i++) {
                bw.write(tokens[i] + DS);
            }
            bw.write(density + DS + visc + DS + TMP + DS + Rtot);
            bw.newLine();
            bw.flush();
        }
        br.close();
        bw.close();
    }

}
