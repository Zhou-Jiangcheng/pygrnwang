import edu.sc.seis.TauP.Arrival;
import edu.sc.seis.TauP.TauP_Time;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;

/** Internal tab-separated protocol; one process handles the complete batch. */
public final class PygrnwangTauPBatch {
    public static void main(String[] args) throws Exception {
        if (args.length < 5) {
            throw new IllegalArgumentException(
                "usage: PygrnwangTauPBatch MODEL SOURCE RECEIVER MODE PHASES...");
        }
        TauP_Time taup = new TauP_Time(args[0]);
        taup.setSourceDepth(Double.parseDouble(args[1]));
        taup.setReceiverDepth(Double.parseDouble(args[2]));
        boolean firstOnly = args[3].equals("first");
        BufferedReader input = new BufferedReader(
            new InputStreamReader(System.in, StandardCharsets.UTF_8));
        PrintWriter output = new PrintWriter(System.out);
        String line;
        int index = 0;
        while ((line = input.readLine()) != null) {
            double distance = Double.parseDouble(line);
            for (int group = 4; group < args.length; group++) {
                taup.setPhaseNames(args[group].isEmpty()
                    ? new String[0] : args[group].split(","));
                taup.calculate(distance);
                int count = taup.getNumArrivals();
                if (count == 0) {
                    output.println(index + "\t" + (group - 4) + "\t\t\tNaN\tNaN");
                }
                for (int arrivalIndex = 0;
                     arrivalIndex < (firstOnly ? Math.min(count, 1) : count);
                     arrivalIndex++) {
                    Arrival arrival = taup.getArrival(arrivalIndex);
                    output.println(index + "\t" + (group - 4) + "\t"
                        + arrival.getName() + "\t" + arrival.getPuristName() + "\t"
                        + Double.toString(arrival.getTime()) + "\t"
                        + Double.toString(arrival.getRayParam()));
                }
            }
            index++;
        }
        output.flush();
    }
}
