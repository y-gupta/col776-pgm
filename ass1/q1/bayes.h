#include <vector>
#include <cstdio>
#include <cassert>
#include <cstring>
#include <iostream>
using namespace std;
class BayesNet{
  public:
    class Node{
    public:
      int id;
      Node *predecessor;
      vector<Node*> parents;
      vector<Node*> children;
      bool observed, descendentObserved, upstream, downstream;
      Node(int _id=0){
        id=_id;
        downstream = false;
        upstream = false;
        predecessor = NULL;
        observed = false;
        descendentObserved = false;
      }
    };
    vector<Node> nodes;
    void load(const char *fname){
      FILE *fs = fopen(fname, "r");
      assert(fs);
      int numNodes = 0;
      assert(1==fscanf(fs,"%d",&numNodes));
      char line[1024];
      assert(fgets(line,1024,fs)); //discard rest of line from first line.

      nodes.resize(0);
      nodes.resize(numNodes);
      for(int i=0;i<numNodes;i++){
        int id, child, bytesRead;
        char *remaining;
        remaining = fgets(line, 1024, fs);
        assert(remaining);
        assert(1==sscanf(remaining,"%d [%n",&id,&bytesRead));
        id--;
        nodes[id].id = id;
        remaining += bytesRead;
        while(1==sscanf(remaining,"%d,%n",&child,&bytesRead)){
          child--;
          nodes[id].children.push_back(&nodes[child]);
          nodes[child].parents.push_back(&nodes[id]);
          remaining += bytesRead;
        }
      }
      fclose(fs);
    }
    void save(const char *fname){
      FILE *ft = fopen(fname, "w");
      assert(ft);
      fprintf(ft,"%d\n",nodes.size());
      for(auto &n: nodes){
        fprintf(ft, "%d [",n.id+1);
        bool first = true;
        for(auto child: n.children){
          if(first)
            fprintf(ft, "%d", child->id+1);
          else
            fprintf(ft, ",%d", child->id+1);
          first = false;
        }
        fprintf(ft, "]\n");
        // first = true;
        // for(auto parent: n.parents){
        //   if(first)
        //     fprintf(ft, "%d", parent->id+1);
        //   else
        //     fprintf(ft, ",%d", parent->id+1);
        //   first = false;
        // }
        // fprintf(ft, "}\n");
      }
      fclose(ft);
    }

    //input query file name, output file name.
    void query(const char *input, const char *output){
      FILE *fs = fopen(input, "r");
      FILE *ft = fopen(output, "w");
      char line[1024];
      while(1){
        int a,b,c, bytesRead;
        vector<int> Z;
        char *remaining;
        remaining = fgets(line, 1024, fs);
        if(!remaining)
          break;
        assert(2==sscanf(remaining,"%d %d [%n",&a,&b,&bytesRead));
        fprintf(ft, "%d %d [",a,b);
        a--;b--;
        remaining += bytesRead;
        while(1==sscanf(remaining,"%d,%n",&c,&bytesRead)){
          if(!Z.empty())
            fprintf(ft, ",");
          fprintf(ft, "%d", c);
          c--;
          Z.push_back(c);
          remaining += bytesRead;
        }
        fprintf(ft, "]\n");
        auto path = activeTrail(b,a,Z);
        if(path.empty()){
          fputs("yes\n",ft);
        }else{
          fprintf(ft, "no [");
          for(int i=0;i<path.size();i++){
            if(i>0)
              fprintf(ft,",");
            fprintf(ft,"%d",1+path[i]);
          }
          fprintf(ft, "]\n");
        }
      }
    }

    //n: numNodes, k=maxChildren
    void generate(int n,int k){
      nodes.resize(0);
      nodes.resize(n);
      for(int i=0;i<n;i++){
        nodes[i].id=i;
        int numNeeded = int(float(rand())/(RAND_MAX)*(k+1));
        for(int j=i+1;j<n;j++){
          if(float(rand())/RAND_MAX <= float(numNeeded)/(n-j)){
            nodes[i].children.push_back(&nodes[j]);
            nodes[j].parents.push_back(&nodes[i]);
            numNeeded--;
          }
        }
      }
    }
    // vector<int> activeTrail1(int a, int b, vector<int> Z, bool isDownstream=true){

    // }
    vector<int> activeTrail(int a, int b, vector<int> Z){
      vector< pair<Node*, bool> > stack;
      for(auto &n: nodes){
        n.upstream = false;
        n.downstream = false;
        n.observed = false;
        n.descendentObserved=false;
        n.predecessor = NULL;
      }
      for(auto i: Z){
        nodes[i].observed = true;
        stack.push_back(make_pair(&nodes[i],false));
      }
      Node *n=NULL;
      while(!stack.empty()){
        n = stack.back().first;
        stack.pop_back();
        n->descendentObserved = true;
        for(auto parent: n->parents){
          if(parent->descendentObserved == false){
            stack.push_back(make_pair(parent,false));
          }
        }
      }
      stack.push_back(make_pair(&nodes[a],false));
      stack.push_back(make_pair(&nodes[a],true));
      nodes[a].downstream = true;
      nodes[a].upstream = true;
      while(!stack.empty()){
        auto node_dir=stack.back();
        auto n = node_dir.first;
        auto isDownstream = node_dir.second;
        // cout<<n->id<<endl;
        stack.pop_back();
        // cout<<n->id<<" ";
        if(n->id == b){
          vector<int> path;
          for(Node *m=n;m != NULL;m = m->predecessor){
            path.push_back(m->id);
          }
          return path;
        }
        for(auto child: n->children){
          if(!n->observed && !child->downstream)
          {
            if(!child->predecessor)
              child->predecessor = n;
            child->downstream = true;
            stack.push_back(make_pair(child,true));
          }
        }
        for(auto parent: n->parents){
          if(isDownstream && n->descendentObserved && !parent->upstream && !parent->observed){
            // v-structure
            if(!parent->predecessor)
              parent->predecessor = n;
            parent->upstream = true;
            stack.push_back(make_pair(parent,false));
          }
          if(!n->observed && !isDownstream && !parent->upstream && !parent->observed)
          {
            if(!parent->predecessor)
              parent->predecessor = n;
            parent->upstream = true;
            stack.push_back(make_pair(parent,false));
          }
        }
      }
      return vector<int>();
    }
};