import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class DNA {

    public static String reverseComplement(String s){
        Map<Character, Character> map = Stream.of(new Character[][] {
                { 'A', 'T' },
                { 'C', 'G' },
                { 'G', 'C' },
                { 'T', 'A' },
                { 'N', 'N' },
        }).collect(Collectors.toMap(data -> data[0], data -> data[1]));
        StringBuilder rev = new StringBuilder();
        for(int i = s.length()-1; i >= 0; i--){
            rev.append(map.get(s.charAt(i)));
        }
        return rev.toString();
    }

    public static String readGenome(String filename){
        StringBuilder contentBuilder = new StringBuilder();
        try (BufferedReader br = new BufferedReader(new FileReader(filename)))
        {
            String sCurrentLine;
            sCurrentLine = br.readLine(); // skip first line
            while ((sCurrentLine = br.readLine()) != null)
            {
                contentBuilder.append(sCurrentLine);
            }
        }
        catch (IOException e)
        {
            e.printStackTrace();
        }
        return contentBuilder.toString();
    }

    public static List<Integer> naive(String p, String t){
        List<Integer> list = new ArrayList<>();
        for(int i = 0 ; i < t.length() - p.length() + 1; i++){
            boolean match = true;
            for(int j = 0 ; j < p.length(); j++){
                if(t.charAt(i+j) != p.charAt(j)){
                    match = false;
                    break;
                }
            }
            if(match){
                list.add(i);
            }
        }
        return list;
    }

    public static List<Integer> naive_with_rc(String p, String t){
        String rc = reverseComplement(p);
        if(rc.equals(p)){
            return naive(p, t);
        }
        List<Integer> list = new ArrayList<>(naive(p,t));
        list.addAll(naive(rc, t));
        return list;
    }


    public static void main(String[] args){
        long start = System.currentTimeMillis();
        String t = readGenome("phix.fa");
        String p = "ATTA";
        List<Integer> occurrences = naive_with_rc(p, t);
        System.out.println("offset of leftmost occurrence: " + Collections.min(occurrences));
        System.out.println("# occurrences: " + occurrences.size());
        long finish = System.currentTimeMillis() - start;
        System.out.println("took: " + finish + " milliseconds");
    }
}
