/*
  merge.cc
  Toshiyuki Gogami, Sep 9, 2019
*/

void merge(){
  TChain*c1 = new TChain("tree");
  c1->Add("tri_Lsingle_111721.root");
  c1->Add("tri_Lsingle_111721_1.root");
  c1->Add("tri_Lsingle_111721_2.root");
  c1->Add("tri_Lsingle_111721_3.root");
  c1->Add("tri_Lsingle_111721_4.root");
  c1->Add("tri_Lsingle_111719.root");
  c1->Add("tri_Lsingle_111719_1.root");
  c1->Add("tri_Lsingle_111719_2.root");
  c1->Merge("test_LHRS.root");
}
