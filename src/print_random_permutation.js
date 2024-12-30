import {random_permutation} from "./src/random.js";
document.open();
document.write("この文章はファイルから表示されています。" + "<br>");
document.write("長さ10のランダムな順列を表示します。" + "<br>");
random_permutation(10).forEach(function(elem) {
    document.write(elem, " ");
});
document.write("<br>");
document.close();