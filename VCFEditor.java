import java.util.*;
import java.io.*;
public class VCFEditor {
static boolean verbose = true;
public static void main(String[] args) throws IOException
{
    String lenFn = args[0], posFn = args[1], vcfFn = args[2];
    String outputFn = args[4];
    String fastaFn = args[3];
    int left = Integer.parseInt(args[5]);
    int right = Integer.parseInt(args[6]);
    Scanner lenInput = new Scanner(new FileInputStream(new File(lenFn)));
    Scanner vcfScanner = new Scanner(new FileInputStream(new File(vcfFn)));
    Scanner posInput = new Scanner(new FileInputStream(new File(posFn)));
    PrintWriter out = new PrintWriter(new File(outputFn));
    PrintWriter logout = new PrintWriter(new File(vcfFn+".log"));
    HashMap<String, Integer> lenMap = new HashMap<String, Integer>();
    String find = "<DEL>";
    String id = "";
    if(lenInput.hasNext())
    {
        id = lenInput.nextLine().substring(1);
    }
    while(lenInput.hasNext())
    {
        String length = lenInput.nextLine();
        if(length.length() > 0 && length.charAt(0) == '>')
        {
            // No line was output for this entry, so save as next id and skip it
            id = length.substring(1);
            continue;
        }
        lenMap.put(id.substring(id.indexOf('_')+1), Integer.parseInt(length));
        while(lenInput.hasNext())
        {
            String cur = lenInput.nextLine();
            if(cur.length() == 0 || cur.charAt(0) != '>') continue;
            id = cur.substring(1);
            break;
        }
    }
    if(verbose) System.out.println("Finished reading " + lenMap.size() + " deletion lengths");
    HashMap<String, Integer> posMap = new HashMap<String, Integer>();
    if(posInput.hasNext())
    {
        id = posInput.nextLine().substring(1);
    }
    while(posInput.hasNext())
    {
        String posString = posInput.nextLine();
        if(posString.length() > 0 && posString.charAt(0) == '>')
        {
            // No line was output for this entry, so save as next id and skip it
            id = posString.substring(1);
            continue;
        }
        else if(posString.length() == 0) continue;
        int p = Integer.parseInt(posString);
        posMap.put(id.substring(id.indexOf('_')+1), p);
        while(posInput.hasNext())
        {
            String cur = posInput.nextLine();
            if(cur.length() == 0 || cur.charAt(0) != '>') continue;
            id = cur.substring(1);
            break;
        }
    }
    int lineNum = 0;
    if(verbose) System.out.println("Finished reading " + posMap.size() + " positions");
    while(vcfScanner.hasNext())
    {
        String line = vcfScanner.nextLine();
        lineNum++;
        if(lineNum%1000 == 0)
        {
            System.out.println("Processed " + lineNum + " lines");
        }
        if(line.charAt(0) == '#')
        {
            out.println(line);
            continue;
        }
        String[] tokens = line.split("\t");
        int index = 0;
        int curPos = Integer.parseInt(tokens[1]);
        String ch = tokens[0];
        if(!line.contains(find))
        {
            out.println(line);
            continue;
        }
        if(!lenMap.containsKey(ch+":"+curPos+""))
        {
            out.println(line);
            logout.println(ch+":"+curPos+"");
            continue;
        }
        Integer deletionLength = lenMap.get(ch+":"+curPos+"");
        Integer svPos = posMap.get(ch+":"+curPos+"");
        String output = line;
        if(deletionLength == 0)
        {
            logout.println("Used sniffles output for " + (ch+":"+curPos+""));
        }
        else
        {
            output = substitute(line, "SVLEN", deletionLength+"");
        }
        StringTokenizer str = new StringTokenizer(output);
        int tokenIdx = 0;
        String seq = "X";
        String newSeq = "";
        if(svPos != null && svPos != -1)
        {
            int actualLeft = Math.min(left, svPos - 1);
            String command = "samtools faidx " + fastaFn + " " + ch + ":" + Math.max(1, (svPos-left+1)) + "-" + (svPos+right+deletionLength);
            if(verbose) System.out.println("deletion: " + command);
            Process child = Runtime.getRuntime().exec(command);
            InputStream seqStream = child.getInputStream();
            Scanner seqInput = new Scanner(seqStream);
            if(!seqInput.hasNext())
            {
                out.println(line);
                logout.println("faidx failed in "+ch+":"+curPos+" "+deletionLength);
                continue;
            }
            seqInput.next();
            if(!seqInput.hasNext())
            {
                out.println(line);
                logout.println("faidx failed in "+ch+":"+curPos+" "+deletionLength);
                continue;
            }
	seq = "";
            while(seqInput.hasNext()) seq += seqInput.next();
            newSeq = seq.substring(0, actualLeft) + seq.substring(actualLeft + deletionLength);
        }
        while(str.hasMoreTokens())
        {
            String cur = str.nextToken();
            if(deletionLength > 0)
            {
                if(svPos != null && svPos != -1 && tokenIdx == 1)
                {
                    cur = svPos+"";
                }
                if(tokenIdx == 3)
                {
                    cur = seq;
                }
                else if(tokenIdx == 4)
                {
                    cur = newSeq;
                }
            }
            out.print(cur);
            if(str.hasMoreTokens()) out.print('\t');
            else out.println();
            tokenIdx++;
        }
    }
    
    out.close();
    logout.close();
}
static String getField(String s, String field)
{
    int n = s.length(), m = field.length();
    for(int i = 0; i+m+1 <= n; i++)
    {
        if(s.substring(i, i+m+1).equals(field+"="))
        {
            String after = s.substring(i+m+1);
            return after.substring(0, after.indexOf(';'));
        }
    }
    return "";
}
static String substitute(String s, String field, String val)
{
    if(val.length() == 0)
    {
        return s;
    }
    int n = s.length(), m = field.length();
    for(int i = 0; i+m+1 <= n; i++)
    {
        if(s.substring(i, i+m+1).equals(field+"="))
        {
            String before = s.substring(0, i+m+1);
            String after = s.substring(i+m+1);
            after = after.substring(after.indexOf(';'));
            return before + val + after;
        }
    }
    return "";
}
}
