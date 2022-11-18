function [BW,maskedImage] = segmentImage(X)
%segmentImage Segment image using auto-generated code from imageSegmenter app
%  [BW,MASKEDIMAGE] = segmentImage(X) segments image X using auto-generated
%  code from the imageSegmenter app. The final segmentation is returned in
%  BW, and a masked image is returned in MASKEDIMAGE.

% Auto-generated by imageSegmenter app on 21-Sep-2021
%----------------------------------------------------


% Normalize input data to range in [0,1].
Xmin = min(X(:));
Xmax = max(X(:));
if isequal(Xmax,Xmin)
    X = 0*X;
else
    X = (X - Xmin) ./ (Xmax - Xmin);
end

% Graph cut
foregroundInd = [3366 5926 8486 10534 12582 14630 17189 19749 21796 22820 25380 26915 28451 29473 31009 31521 33056 34080 35614 37150 39710 40734 42270 44318 45342 46878 48414 50974 51998 53534 55070 56606 57630 59166 59678 60702 62239 63263 63775 64799 65312 66336 66848 68898 70434 71971 72483 73509 75045 76070 78118 78631 80168 82730 83754 87853 88365 89902 90414 91438 92976 93488 95537 96050 97074 98099 99123 100661 101686 102710 103224 103736 104760 106297 106811 107323 108348 108860 109373 110398 110910 111936 112448 112961 113473 114499 115011 116036 116548 117572 117574 118086 118598 118599 ];
backgroundInd = [13 527 819 1044 1331 1696 1697 1699 1702 1705 1709 1712 1715 1718 1720 1723 1724 1727 1728 1730 1732 1735 1737 1739 1742 1745 1916 1919 1922 1925 1929 1931 1934 2073 2195 2200 2203 2205 2258 2261 2264 2266 2270 2355 2419 2421 2423 2426 2450 2453 2457 2463 2467 2472 2479 2785 2789 2792 2995 3000 3095 3096 3097 3107 3111 3207 3213 3308 3314 3437 3440 3520 3523 3528 3531 3703 3708 3712 3716 3945 3948 4139 4141 4146 4149 4151 4157 4160 4165 4169 4171 4176 4180 4182 4187 4191 4195 4199 4202 4204 4207 4212 4340 4343 4345 4347 4348 4350 4403 4560 4563 4566 4568 4569 4572 4575 4577 4631 4863 4865 4868 4869 4870 4871 4874 4965 4967 5089 5091 5092 5094 5095 5097 5098 5100 5435 5476 5477 5900 5901 5903 5904 5948 5986 5988 6126 6127 6129 6130 6167 6290 6416 6451 6471 6498 6644 6645 6647 6822 6928 7159 7161 7162 7191 7440 7484 7494 7522 7674 7676 8004 8033 8215 8499 8700 8823 8976 9022 9027 9057 9213 9537 9752 10035 10264 10512 10558 10560 10749 10777 11072 11103 11288 11571 11582 12048 12095 12127 12824 12985 13107 13117 13119 13409 13628 13631 13849 14096 14361 14655 14687 14847 15161 15385 15632 16222 16408 16696 16921 17167 17203 17246 17718 17842 18191 18456 18703 18781 18859 18869 18968 19253 19763 19967 20239 20317 20787 20829 20896 21528 21774 21811 22091 22365 22554 22834 22968 23064 23389 23822 24413 24575 24881 24924 25358 25395 25801 26005 26460 26652 26893 27441 27484 27670 28429 28467 28508 28671 28977 29624 29965 30003 30044 30093 30513 31426 31916 32049 32092 32285 32525 33301 33434 33585 34008 34061 34113 34302 34624 34625 34652 34693 35595 35633 35648 36160 36920 37002 37184 37681 37724 38155 38208 38277 39217 39445 39453 39659 39691 39813 39931 40241 40256 40284 41227 41352 41597 41777 43313 43356 43408 43787 43840 43933 43945 43954 44538 44849 44978 45323 45496 45916 46000 46613 46897 47647 47829 47883 47923 48375 48960 48988 49656 49932 50289 50483 52572 52653 52783 53165 53268 53517 53555 53751 54592 55132 55565 56115 56822 57140 57376 57545 57546 57613 57692 58688 59083 59188 59380 59663 59924 60726 61182 61427 61788 62223 62262 62389 62774 62784 62888 63078 63251 63759 64275 64311 64530 64787 65011 65811 65847 65856 65884 66319 67347 67384 67856 68082 68113 68371 68408 68640 68883 68928 68956 69392 69946 70345 70362 70421 70640 71696 71726 71952 71994 72887 72981 73200 73488 73744 73911 74005 74043 74076 75247 75792 76049 76054 76092 76895 77328 77590 77807 78097 78369 78557 78614 78655 78685 78757 79343 79376 79580 79640 80192 80604 80657 81176 81728 81903 82448 82712 83266 83294 84181 84223 84249 84290 84398 84753 84975 85273 85827 86228 86810 87058 87216 87363 87826 88097 88112 88347 88416 88901 89285 89583 89885 89925 91422 91924 92486 92766 92960 93025 93169 93387 93405 94023 94406 94432 94497 95265 95445 95508 95521 95561 95653 95730 96586 96963 97059 97099 97121 98290 98596 98638 99093 99345 99503 99621 101054 101157 101158 101200 101219 101607 101653 101876 101905 102226 102433 102696 103599 103605 103701 103720 103989 104233 104276 104703 105259 105489 105771 106485 106775 106796 106837 106852 106917 107040 107820 107821 107862 108227 108333 108963 109152 109250 109359 109398 109557 109924 110359 110623 110936 111408 112349 112473 112484 112670 112920 113143 113163 113459 114718 114968 115547 116215 116532 117084 117092 117160 118040 118070 118108 118109 118301 118775 119645 119652 120558 120631 121181 121371 121626 121849 122212 122682 122717 122794 122939 123399 123742 123930 124731 124766 124772 124921 125111 125243 125688 125723 126268 126894 127002 127591 127705 127993 128356 128701 129440 129562 129820 130997 131578 131938 132634 134075 134429 134498 134682 135674 136371 136545 137242 137501 137793 138176 138593 138961 139260 139960 140575 140826 141665 142003 143136 143356 143711 144051 144410 145644 145698 145759 145827 145861 146722 147572 147806 148172 148259 148476 148771 149530 150309 150602 151282 151333 151389 152060 152358 152777 153895 154138 154973 155432 156873 156970 157482 157610 157637 157693 157897 158507 158804 159068 159435 159770 161069 161497 161789 162138 162606 162947 163029 163632 164274 165169 165209 165373 165404 166493 166706 167421 168243 168279 168375 168672 168896 168957 169269 169981 170326 170806 170836 170941 171318 171517 171625 172224 172371 172480 172573 172856 173053 173394 173881 173905 174395 174414 174589 174715 174908 174909 174910 174924 175274 175424 175764 175947 176269 176370 176372 176448 176449 176456 176637 176669 176961 177479 177985 177990 178106 178499 178500 179197 180994 181249 181791 182781 184082 184829 185375 185632 186160 187047 187309 187709 187936 188413 189769 191008 191316 191485 192350 192929 193056 193896 194557 194744 195440 195616 195826 195960 195967 195978 195989 197117 197125 197664 199165 199883 199914 200189 200924 201248 201725 203261 203773 204832 205309 205831 206845 207904 208381 209405 209917 210464 210941 212477 212512 214013 214025 214560 215037 215072 215549 216096 216573 216809 216815 216823 216829 216837 216843 216848 216853 216857 216861 216864 217379 217384 217632 217806 217812 217819 217826 218108 218309 218928 219146 219168 219309 219320 219448 219644 219970 220322 220704 220824 221003 221180 221196 221525 221728 222204 222221 222240 222348 222558 222716 222736 223078 223264 223740 223763 223871 224252 224276 224620 224626 224792 224795 224883 225146 225276 225311 225314 225662 225717 225731 225896 226219 226299 226300 226342 226345 226692 226696 226700 226704 226711 226717 226723 226767 226788 226800 226808 226858 227412 227421 227888 227894 227900 227908 227915 ];
L = superpixels(X,841);
BW = lazysnapping(X,L,foregroundInd,backgroundInd);

% Local graph cut
xPos = [8.9000 27.1857 51.3229 82.7743 128.1229 155.1857 204.9229 237.8371 253.1971 243.6886 212.2371 168.3514 150.7971 112.0314 76.1914 21.3343 8.9000 5.9743 ];
yPos = [321.9629 304.4086 300.0200 295.6314 300.7514 313.1857 325.6200 336.5914 327.0829 315.3800 305.1400 287.5857 282.4657 279.5400 276.6143 283.1971 288.3171 313.1857 ];
m = size(BW, 1);
n = size(BW, 2);
ROI = poly2mask(xPos,yPos,m,n);
foregroundInd = [];
backgroundInd = [];
L = superpixels(X,1257);
BW = BW | grabcut(X,L,ROI,foregroundInd,backgroundInd);

% Active contour
iterations = 100;
BW = activecontour(X, BW, iterations, 'Chan-Vese');

% Create masked image.
maskedImage = X;
maskedImage(~BW) = 0;
end

