import {random_permutation} from "./random.js";

document.write("長さ10のランダムな順列を表示します。" + "<br>");

random_permutation(10).forEach(function(elem) {
    document.write(elem, " ");
});
document.write("<br>");